#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/repdiag.h>

#include <iostream>
#include <cstdlib>
#include <cmath>

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkMassProperties.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCell.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkIdList.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkImplicitPolyDataDistance.h>
Eigen::MatrixXd V,U;
Eigen::MatrixXi F;
Eigen::SparseMatrix<double> L;

int iter;
int iter2;
char temp[128];
char temp2[128];

int main(int argc, char *argv[]) {
    using namespace Eigen;
    using namespace std;

    if (argc != 2) {
        cerr << "Usage " << argv[0] <<": inputMesh" << endl;
        return -1;
    }
    igl::readOFF(argv[1], V, F);
    // construct polys
    vtkSmartPointer<vtkCellArray> polys =
        vtkSmartPointer<vtkCellArray>::New();
    for(int i = 0; i < F.rows(); ++i) {
        vtkSmartPointer<vtkIdList> ids =
            vtkSmartPointer<vtkIdList>::New();
        ids->InsertNextId(F(i,0));
        ids->InsertNextId(F(i,1));
        ids->InsertNextId(F(i,2));
        polys->InsertNextCell(ids);
    }

    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
        vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    smoother->SetNumberOfIterations(50);
    smoother->BoundarySmoothingOff();
    smoother->FeatureEdgeSmoothingOff();
    smoother->SetPassBand(0.001);
    smoother->NonManifoldSmoothingOn();
    smoother->NormalizeCoordinatesOn();

    vtkSmartPointer<vtkMassProperties> mass =
        vtkSmartPointer<vtkMassProperties>::New();
    //    writer->SetFileName("temp.vtk");
    //    writer->SetInputData(polydata);
    //    writer->Update();
    igl::cotmatrix(V,F,L);
    U = V;
    double q = 1.0; // for now 0.2
    while(q > 0.1) {
        // construct points
        cout << "iteration " << iter << endl;
        vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();
        for(int i = 0; i < U.rows(); ++i) {
            points->InsertNextPoint(U(i,0), U(i,1), U(i,2));
        }
        vtkSmartPointer<vtkPolyData> polydata =
            vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);
        polydata->SetPolys(polys);
        // smooth it
        smoother->SetInputData(polydata);
        smoother->Update();
        vtkSmartPointer<vtkPolyData> polydata_smooth = smoother->GetOutput();
        mass->SetInputData(polydata_smooth);
        mass->Update();
        double current_volume = mass->GetVolume();
        for(int i = 0; i < U.rows(); ++i) {
            double p[3];
            polydata_smooth->GetPoint(i,p);
            U(i,0) = p[0];
            U(i,1) = p[1];
            U(i,2) = p[2];
        }
        SparseMatrix<double> M;
        igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
        // Solve (M-delta*L) U = M*U
        const auto & S = (M - 1e-4*L);
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
        assert(solver.info() == Eigen::Success);
        U = solver.solve(M*U).eval();
        // Compute centroid and subtract (also important for numerics)
        VectorXd dblA;
        igl::doublearea(U,F,dblA);
        double area = 0.5*dblA.sum();
        MatrixXd BC;
        igl::barycenter(U,F,BC);
        RowVector3d centroid(0,0,0);
        for(int i = 0;i<BC.rows();i++)
        {
            centroid += 0.5*dblA(i)/area*BC.row(i);
        }
        U.rowwise() -= centroid;
        // Normalize to unit surface area (important for numerics)
        U.array() /= sqrt(area);
        sprintf(temp,"temp_%04d.off", iter++);
        // get the ellipsoid from the second moment matrix
        // U is n by 3
        MatrixXd cog = U.colwise().mean();
        MatrixXd U_centered = U - cog.replicate(U.rows(), 1);
        MatrixXd U_transposed = U_centered.transpose();
        MatrixXd second_moment = U_transposed * U_centered;

        SelfAdjointEigenSolver<MatrixXd> es(second_moment);
        VectorXd radii = es.eigenvalues();
        cout << radii << endl;
        MatrixXd rotation = es.eigenvectors(); // 3 by 3 rotation matrix

        double ellipsoid_volume = 4 / 3.0 * M_PI * sqrt(radii(0)) * sqrt(radii(1)) * sqrt(radii(2));
        double volume_factor = pow(current_volume / ellipsoid_volume, 1.0 / 3.0);
        // obtain the best fitting ellipsoid from the second moment matrix
        vtkSmartPointer<vtkParametricEllipsoid> ellipsoid =
            vtkSmartPointer<vtkParametricEllipsoid>::New();
        ellipsoid->SetXRadius(sqrt(radii(0))*volume_factor);
        ellipsoid->SetYRadius(sqrt(radii(1))*volume_factor);
        ellipsoid->SetZRadius(sqrt(radii(2))*volume_factor);

        vtkSmartPointer<vtkParametricFunctionSource> parametric_function =
            vtkSmartPointer<vtkParametricFunctionSource>::New();
        parametric_function->SetParametricFunction(ellipsoid);
        parametric_function->Update();

        vtkSmartPointer<vtkPolyData> ellipsoid_polydata = parametric_function->GetOutput();
        int num_points = ellipsoid_polydata->GetNumberOfPoints();
        // n by 3
        MatrixXd ellipsoid_points_matrix(num_points, 3);
        for(int i = 0; i < num_points; ++i) {
            double p[3];
            ellipsoid_polydata->GetPoint(i, p);
            ellipsoid_points_matrix(i,0) = p[0];
            ellipsoid_points_matrix(i,1) = p[1];
            ellipsoid_points_matrix(i,2) = p[2];
            if(abs(p[0]) > sqrt(radii[0])*volume_factor || abs(p[1]) > sqrt(radii[1])*volume_factor || abs(p[2]) > sqrt(radii[2]) * volume_factor ) {
                cerr << "Wrong ellipsoid point" << endl << endl;
                cout << "(" << p[0] << "," << p[1] << "," << p[2] << ")" << endl;
                cout << ellipsoid->GetXRadius() << endl;
                cout << ellipsoid->GetYRadius() << endl;
                cout << ellipsoid->GetZRadius() << endl;
            }
        }
        // 3 by n
        ellipsoid_points_matrix.transposeInPlace();
        MatrixXd rotated_points = rotation * ellipsoid_points_matrix;
        // n by 3
        rotated_points.transposeInPlace();
        MatrixXd translated_rotated_points = rotated_points + cog.replicate(num_points, 1);

        vtkSmartPointer<vtkPoints> ellipsoid_points =
            vtkSmartPointer<vtkPoints>::New();
        for(int i = 0; i < num_points; ++i) {
            ellipsoid_points->InsertNextPoint(translated_rotated_points(i,0), translated_rotated_points(i,1), translated_rotated_points(i,2));
            //cout << "(" << translated_rotated_points(i,0) << "," << translated_rotated_points(i,1) << "," << translated_rotated_points(i,2) << ")" << endl;
        }
        ellipsoid_polydata->SetPoints(ellipsoid_points);
        ellipsoid_polydata->Modified();

        vtkSmartPointer<vtkImplicitPolyDataDistance> implicit_distance_filter = 
            vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
        implicit_distance_filter->SetInput(ellipsoid_polydata);
        vector<float> distance_vector;
        for(int i = 0; i < U.rows(); ++i) {
            double p[3];
            p[0] = U(i,0);
            p[1] = U(i,1);
            p[2] = U(i,2);
            float current_distance = abs(implicit_distance_filter->EvaluateFunction(p));
            distance_vector.push_back(current_distance);
        }
        sort(distance_vector.begin(), distance_vector.end());
        string fname = temp;
        igl::writeOFF(fname, U, F);
        int quantile_idx = static_cast<int>(U.rows() * 0.95);
        cout << distance_vector[0] << endl;
        cout << distance_vector[U.rows()-1] << endl;

        q = distance_vector[quantile_idx-1];
        sprintf(temp2,"ell_%04d.vtk", iter2++);
        string fname2 = temp2;
        writer->SetFileName(temp2);
        writer->SetInputData(ellipsoid_polydata);
        writer->Update();
        cout << endl << endl;
    }
    return 0;
}
