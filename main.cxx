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
#include <igl/readSTL.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/repdiag.h>
#include <igl/opengl/glfw/Viewer.h>


#include <iostream>
#include <fstream>
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
#include <vtkSmoothPolyDataFilter.h>
#include <vtkIdList.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkCleanPolyData.h>

igl::opengl::glfw::Viewer viewer;

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[]) {

    if (argc != 5) {
        cerr << "Usage " << argv[0] <<": inputMesh dt smoothAmount maxIter" << endl;
        return -1;
    }


    MatrixXd temp_V, V, U;
    MatrixXi temp_F, F, SVI, SVJ;
    MatrixXd N;
    SparseMatrix<double> L;
    MatrixXd rotation;
    VectorXd radii;

    double dt = atof(argv[2]);
    double smooth_amount = atof(argv[3]);
    ofstream radii_file;
    radii_file.open("radii.txt");
    int iter = 0;
    int max_iter = atoi(argv[4]);
    printf("Opening %s\n",argv[1]);
    igl::readSTL(argv[1], temp_V, temp_F, N);
    igl::remove_duplicate_vertices(temp_V,temp_F,0,V,SVI,SVJ,F);
    printf("Read STL file\n");
    igl::cotmatrix(V,F,L);
    // cout << "calculated  L "<< endl << L << endl;
    U = V;

    double q = 1.0; // for now 0.2

    const auto &key_down = [&q, &U, &V, &F, &dt, &L, &iter, &smooth_amount,&rotation, &radii, &max_iter, &radii_file](igl::opengl::glfw::Viewer &viewer,unsigned char key,int mod)->bool
    {
      switch(key)
      {
        case 'r':
        case 'R':
          U = V;
          break;
        case ' ':
        {
          if(q > 0.0) {
              char temp[128];
              // compute mean curvature flow
              SparseMatrix<double> M;
              igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
              // Solve (M-delta*L) U = M*U
              const auto & S = (M - dt*L);
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

              vtkSmartPointer<vtkPoints> points =
                  vtkSmartPointer<vtkPoints>::New();
              for(int i = 0; i < U.rows(); ++i) {
                  points->InsertNextPoint(U(i,0), U(i,1), U(i,2));
              }

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

              vtkSmartPointer<vtkPolyData> polydata =
                  vtkSmartPointer<vtkPolyData>::New();
              polydata->SetPolys(polys);
              polydata->SetPoints(points);

              // polydata writer
              vtkSmartPointer<vtkPolyDataWriter> writer =
                  vtkSmartPointer<vtkPolyDataWriter>::New();

              // smoother
              vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
                  vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
              smoother->SetNumberOfIterations(20);
              smoother->BoundarySmoothingOff();
              smoother->FeatureEdgeSmoothingOff();
              smoother->SetPassBand(smooth_amount);
              smoother->NonManifoldSmoothingOn();
              smoother->NormalizeCoordinatesOn();

              // mass filter
              vtkSmartPointer<vtkMassProperties> mass =
                  vtkSmartPointer<vtkMassProperties>::New();
              mass->SetInputData(polydata);
              mass->Update();

              // smooth polydata
              smoother->SetInputData(polydata);
              smoother->Update();
              vtkSmartPointer<vtkPolyData> polydata_smooth = smoother->GetOutput();
              //mass->SetInputData(polydata_smooth);
              //mass->Update();
              double current_volume = mass->GetVolume();

              // set U to smoothed points
              for(int i = 0; i < U.rows(); ++i) {
                  double p[3];
                  polydata->GetPoint(i,p);
                  U(i,0) = p[0];
                  U(i,1) = p[1];
                  U(i,2) = p[2];
              }

              vtkSmartPointer<vtkCleanPolyData> cleaner =
                  vtkSmartPointer<vtkCleanPolyData>::New();
              cleaner->SetTolerance(0.05);
              cleaner->SetInputData(polydata_smooth);
              cleaner->Update();
              vtkSmartPointer<vtkPolyData> cleaned_polydata = cleaner->GetOutput();
              MatrixXd U_temp(cleaned_polydata->GetNumberOfPoints(),3);
              // set U to smoothed points
              for(int i = 0; i < U_temp.rows(); ++i) {
                  double p[3];
                  cleaned_polydata->GetPoint(i,p);
                  U_temp(i,0) = p[0];
                  U_temp(i,1) = p[1];
                  U_temp(i,2) = p[2];
              }
              // filter out clustered points
              MatrixXd cog = U_temp.colwise().mean();
              MatrixXd U_centered = U_temp - cog.replicate(U_temp.rows(), 1); // N by 3
              MatrixXd U_transposed = U_centered.transpose();
              Matrix3d U_second_moment = U_transposed * U_centered;

              SelfAdjointEigenSolver<MatrixXd> es(U_second_moment);
              rotation = es.eigenvectors(); // 3 by 3 rotation matrix
              radii = es.eigenvalues();
              radii(0) = sqrt(radii(0));
              radii(1) = sqrt(radii(1));
              radii(2) = sqrt(radii(2));

              double ellipsoid_volume = 4 / 3.0 * M_PI * radii(0) * radii(1) * radii(2);
              double volume_factor = pow(current_volume / ellipsoid_volume, 1.0 / 3.0);
              radii(0) *= volume_factor;
              radii(1) *= volume_factor;
              radii(2) *= volume_factor;
              // obtain the best fitting ellipsoid from the second moment matrix
              vtkSmartPointer<vtkParametricEllipsoid> ellipsoid =
                  vtkSmartPointer<vtkParametricEllipsoid>::New();
              ellipsoid->SetXRadius(radii(0));
              ellipsoid->SetYRadius(radii(1));
              ellipsoid->SetZRadius(radii(2));

              vtkSmartPointer<vtkParametricFunctionSource> parametric_function =
                  vtkSmartPointer<vtkParametricFunctionSource>::New();
              parametric_function->SetParametricFunction(ellipsoid);
              parametric_function->SetUResolution(30);
              parametric_function->SetVResolution(30);
              parametric_function->Update();
              vtkSmartPointer<vtkPolyData> ellipsoid_polydata = parametric_function->GetOutput();

              // Get ellipsoid points into the matrix
              MatrixXd ellipsoid_points_matrix(ellipsoid_polydata->GetNumberOfPoints(), 3);
              for(int i = 0; i < ellipsoid_polydata->GetNumberOfPoints(); ++i) {
                  double p[3];
                  ellipsoid_polydata->GetPoint(i,p);
                  ellipsoid_points_matrix(i,0) = p[0];
                  ellipsoid_points_matrix(i,1) = p[1];
                  ellipsoid_points_matrix(i,2) = p[2];
              }
              // rotate the points
              MatrixXd rotated_ellipsoid_points = rotation * (ellipsoid_points_matrix.transpose());
              rotated_ellipsoid_points.transposeInPlace(); // n x 3
              // translate the points
              MatrixXd translated_points = rotated_ellipsoid_points + cog.replicate(rotated_ellipsoid_points.rows(),1);

              // convert eigen matrix to vtk polydata
              vtkSmartPointer<vtkPolyData> best_fitting_ellipsoid_polydata =
                  vtkSmartPointer<vtkPolyData>::New();
              vtkSmartPointer<vtkPoints> best_fitting_ellipsoid_points =
                  vtkSmartPointer<vtkPoints>::New();
              for(int i = 0; i < translated_points.rows(); ++i) {
                  double p[3] = {translated_points(i,0), translated_points(i,1), translated_points(i,2)};
                  best_fitting_ellipsoid_points->InsertNextPoint(p);
              }
              best_fitting_ellipsoid_polydata->SetPoints(best_fitting_ellipsoid_points);
              best_fitting_ellipsoid_polydata->SetPolys(ellipsoid_polydata->GetPolys());
              best_fitting_ellipsoid_polydata->Modified();

              vtkSmartPointer<vtkImplicitPolyDataDistance> implicit_distance_filter =
                  vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
              implicit_distance_filter->SetInput(best_fitting_ellipsoid_polydata);
              vector<double> distance_vector;

              for(int i = 0; i < U.rows(); ++i) {
                  double p[3];
                  polydata->GetPoint(i,p);
                  double current_distance = abs(implicit_distance_filter->EvaluateFunction(p));
                  distance_vector.push_back(current_distance);
              }

              sort(distance_vector.begin(), distance_vector.end());
              string off_filename = "temp_off/temp_" + prefix;
              off_filename += ".off";
              igl::writeOFF(off_filename, U, F);
              int quantile_idx = static_cast<int>(U.rows() * 0.95);

              q = distance_vector[quantile_idx];
              cout << "iter " << iter << ": " << q << endl;
              string vtk_filename = "ell/ell_" + prefix;
              vtk_filename += ".vtk";
              writer->SetFileName(vtk_filename.c_str());
              writer->SetInputData(best_fitting_ellipsoid_polydata);
              writer->Update();

              vtk_filename = "temp_vtk/temp_" + prefix;
              vtk_filename += ".vtk";
              writer->SetFileName(vtk_filename.c_str());
              writer->SetInputData(polydata_smooth);
              writer->Update();
              string cleaned_vtk_filename = "cleaned/cleaned_" + prefix;
              cleaned_vtk_filename+=".vtk";
              writer->SetFileName(cleaned_vtk_filename.c_str());
              writer->SetInputData(cleaned_polydata);
              writer->Update();
              char buffer[128];
              sprintf(buffer, "%f, %f, %f\n", radii(0)*volume_factor, radii(1) * volume_factor, radii(2) * volume_factor);
              radii_file << buffer;
              if(iter > max_iter) {
                  cerr << "Failed to converge\n" << endl;
                  break;
              }
          }

          break;
        }
        default:
          return false;
      }
      // Send new positions, update normals, recenter
      viewer.data().set_vertices(U);
      viewer.data().compute_normals();
      viewer.core.align_camera_center(U,F);
      return true;
    };

    // Use original normals as pseudo-colors
    igl::per_vertex_normals(V,F,N);
    MatrixXd C = N.rowwise().normalized().array()*0.5+0.5;

    U=V;
    viewer.data().set_mesh(U, F);
    viewer.data().set_colors(C);
    viewer.callback_key_down = key_down;


    cout<<"Press [space] to smooth."<<endl;;
    cout<<"Press [r] to reset."<<endl;;
    viewer.launch();


    cout << radii << endl;
    cout << endl;
    cout << rotation << endl;
    radii_file.close();
}
