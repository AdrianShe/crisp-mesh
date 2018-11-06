#include "marching_cubes_offset.h"
#include "random_points_offset.h"
// #include "poisson_surface_reconstruction.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>


Eigen::MatrixXd V, V_off, V_p, P, N;
Eigen::MatrixXi F, F_off, F_p;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
     if (key == '1') {
        std::cout << "original" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V, F);
  viewer.data().set_points(P, Eigen::RowVector3d(1,1,1));
	viewer.core.align_camera_center(V, F);
     }
     else if (key == '2') {
        std::cout << "offset" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V_off, F_off);
	viewer.core.align_camera_center(V_off, F_off);
     }
     else if (key == '3') {
	std::cout << "offset using random points" << std::endl;
	viewer.data().clear();
 	// viewer.data().set_mesh(V_p, F_p);
	viewer.core.align_camera_center(V_off, F_off);
	// viewer.data().set_points(P, Eigen::RowVector3d(1,1,1));
     }
     return false;
}

int main(int argc, char *argv[])
{
  // Load mesh with desired distance and resolution
  // argv[1] mesh
  // argv[2] distance
  // argv[3] resolution

  igl::read_triangle_mesh(argv[1], V, F);
  double sigma = atof(argv[2]);
  int res = atol(argv[3]);
  std::cout << "distance" << sigma << "resolution" << res << std::endl;
  marching_cubes_offset(V, F, sigma, res, V_off, F_off);
  random_points_offset(V, F, sigma, V.rows() / 20, P, N);
  // poisson_surface_reconstruction(P, N, V_p, F_p);

  // Render the original mesh and the offset mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);
 // viewer.data().set_points(P, Eigen::RowVector3d(1,1,1));
  viewer.launch();
  return EXIT_SUCCESS;
}

