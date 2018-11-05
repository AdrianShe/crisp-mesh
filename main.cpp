#include "marching_cubes_offset.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>


Eigen::MatrixXd V, V_off;
Eigen::MatrixXi F, F_off;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
     if (key == '1') {
        std::cout << "original" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V, F);
	viewer.core.align_camera_center(V, F);
     }
     else if (key == '2') {
        std::cout << "offset" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V_off, F_off);
	viewer.core.align_camera_center(V_off, F_off);
     }
     return false;
}

int main(int argc, char *argv[])
{
  // Load mesh and desired distance
  // argv[1] mesh
  // argv[2] distance
  igl::read_triangle_mesh(argv[1], V, F);
  double sigma = atof(argv[2]);
  int res = atol(argv[3]);
  std::cout << sigma << res << std::endl;
  marching_cubes_offset(V, F, sigma, res, V_off, F_off);
  std::cout << "done" << V_off.cols() << F_off.cols() << std::endl;

  // Render the original mesh and the offset mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);
  viewer.launch();
  return EXIT_SUCCESS;
}

