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
char temp[128];

int main(int argc, char *argv[]) {
    using namespace Eigen;
    using namespace std;

    if (argc != 2) {
        cerr << "Usage " << argv[0] <<": inputMesh" << endl;
        return -1;
    }
    igl::readOFF(argv[1], V, F);

	// construct polys from F
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
	// declare polydata
    vtkSmartPointer<vtkPolyData> polydata =
        vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPolys(polys);

	// polydata writer
    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
	// smoother
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
        vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    smoother->SetNumberOfIterations(20);
    smoother->BoundarySmoothingOff();
    smoother->FeatureEdgeSmoothingOff();
    smoother->SetPassBand(0.01);
    smoother->NonManifoldSmoothingOn();
    smoother->NormalizeCoordinatesOn();

	// mass filter
    vtkSmartPointer<vtkMassProperties> mass =
        vtkSmartPointer<vtkMassProperties>::New();

    igl::cotmatrix(V,F,L);

	MatrixXd cog = V.colwise().mean();
	MatrixXd V_centered = V - cog.replicate(V.rows(), 1); // N by 3
	MatrixXd V_transposed = V_centered.transpose();
	Matrix3d V_second_moment = V_transposed * V_centered;

	SelfAdjointEigenSolver<MatrixXd> es(V_second_moment);
	Matrix3d rotation = es.eigenvectors(); // 3 by 3 rotation matrix

	V_centered.transposeInPlace(); // 3 by N
	rotation.transposeInPlace(); // rotation transposed
	MatrixXd V_temp =  rotation * V_centered; // rotated U
	V_temp.transposeInPlace();

	U = V_temp;
	vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();
	for(int i = 0; i < U.rows(); ++i) {
		points->InsertNextPoint(U(i,0), U(i,1), U(i,2));
	}

	polydata->SetPoints(points);
	polydata->Modified();
	mass->SetInputData(polydata);
	mass->Update();
	double current_volume = mass->GetVolume();
	double scale_factor = pow(current_volume, 1.0 / 3.0);
	vtkSmartPointer<vtkPoints> points_volume_normalized =
			vtkSmartPointer<vtkPoints>::New();
	for(int i = 0; i < U.rows(); ++i) {
		U(i,0) *= scale_factor;
		U(i,1) *= scale_factor;
		U(i,2) *= scale_factor;
		points_volume_normalized->InsertNextPoint(U(i,0), U(i,1), U(i,2));
	}
	polydata->SetPoints(points_volume_normalized);
	polydata->Modified();

	double q = 1.0; // for now 0.2
    while(q > 0.1) {
        // compute mean curvature flow
        SparseMatrix<double> M;
        igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
        // Solve (M-delta*L) U = M*U
        const auto & S = (M - 0.001*L);
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
        sprintf(temp,"%04d", iter++);
	    string prefix = temp;

	    // get the flowed points from U matrix
        vtkSmartPointer<vtkPoints> flowed_points =
            vtkSmartPointer<vtkPoints>::New();
        for(int i = 0; i < U.rows(); ++i) {
	        flowed_points->InsertNextPoint(U(i,0), U(i,1), U(i,2));
        }
        polydata->SetPoints(flowed_points);
        polydata->Modified();

	    // smooth polydata
        smoother->SetInputData(polydata);
        smoother->Update();
        vtkSmartPointer<vtkPolyData> polydata_smooth = smoother->GetOutput();
        mass->SetInputData(polydata_smooth);
//	    mass->SetInputData(polydata);
        mass->Update();

	    // set U to smoothed points
        for(int i = 0; i < U.rows(); ++i) {
            double p[3];
            polydata_smooth->GetPoint(i,p);
            U(i,0) = p[0];
            U(i,1) = p[1];
            U(i,2) = p[2];
        }

        MatrixXd U_transposed = U.transpose();
        MatrixXd U_second_moment = U_transposed * U;

        es.compute(U_second_moment);
        VectorXd radii = es.eigenvalues();
	    // thought I could just call sqrt on the vector...
	    radii(0) = sqrt(radii(0));
	    radii(1) = sqrt(radii(1));
	    radii(2) = sqrt(radii(2));

        double ellipsoid_volume = 4 / 3.0 * M_PI * radii(0) * radii(1) * radii(2);
        double volume_factor = pow(1 / ellipsoid_volume, 1.0 / 3.0); // unit volume
        // obtain the best fitting ellipsoid from the second moment matrix
        vtkSmartPointer<vtkParametricEllipsoid> ellipsoid =
            vtkSmartPointer<vtkParametricEllipsoid>::New();
        ellipsoid->SetXRadius(radii(0)*volume_factor);
        ellipsoid->SetYRadius(radii(1)*volume_factor);
        ellipsoid->SetZRadius(radii(2)*volume_factor);

        vtkSmartPointer<vtkParametricFunctionSource> parametric_function =
            vtkSmartPointer<vtkParametricFunctionSource>::New();
        parametric_function->SetParametricFunction(ellipsoid);
        parametric_function->Update();

        vtkSmartPointer<vtkPolyData> ellipsoid_polydata = parametric_function->GetOutput();

        vtkSmartPointer<vtkImplicitPolyDataDistance> implicit_distance_filter = 
            vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
        implicit_distance_filter->SetInput(ellipsoid_polydata);
        vector<double> distance_vector;
        for(int i = 0; i < U.rows(); ++i) {
            double p[3];
	        polydata->GetPoints()->GetPoint(i,p);
            double current_distance = abs(implicit_distance_filter->EvaluateFunction(p));
            distance_vector.push_back(current_distance);
        }
        sort(distance_vector.begin(), distance_vector.end());
        string off_filename = "temp_" + prefix;
	    off_filename += ".off";
        igl::writeOFF(off_filename, U, F);
        int quantile_idx = static_cast<int>(U.rows() * 0.95);

        q = distance_vector[quantile_idx-1];
	    cout << "iter " << iter << ": " << q << endl;
        string vtk_filename = "ell_" + prefix;
	    vtk_filename += ".vtk";
        writer->SetFileName(vtk_filename.c_str());
        writer->SetInputData(ellipsoid_polydata);
        writer->Update();

	    vtk_filename = "temp_" + prefix;
	    vtk_filename += ".vtk";
	    writer->SetFileName(vtk_filename.c_str());
	    writer->SetInputData(polydata_smooth);
	    writer->Update();
    }
    return 0;
}
