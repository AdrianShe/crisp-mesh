#include "marching_cubes_offset.h"
#include "random_points_offset.h"
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/offset_surface.h>
#include <igl/opengl/glfw/Viewer.h>
#include "igl/jet.h"
#include <Eigen/Core>
#include <string>
#include <iostream>
#include "validate.h"


Eigen::MatrixXd V, V_mc, V_mcr, V_igl, closest_points, C;
Eigen::MatrixXi F, F_mc, F_mcr, F_igl;
Eigen::VectorXd int_dist_o, int_dist_n;
double sigma;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
     if (key == '1') {
        std::cout << "original" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V, F);
	viewer.core.align_camera_center(V, F);
     }
     else if (key == '2') {
        std::cout << "offset using marching cubes" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V_mc, F_mc);
	igl::jet(int_dist_o, 0.5 * sigma, 1.5 * sigma, C);
	viewer.data().set_colors(C);
	// viewer.core.align_camera_center(V_mc, F_mc);
     }
     else if (key == '3') {
	std::cout << "offset using marching cubes and root finding" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V_mcr, F_mcr);
	igl::jet(int_dist_n, 0.5 * sigma, 1.5 * sigma, C);
	viewer.data().set_colors(C);
	// viewer.core.align_camera_center(V_mcr, F_mcr);
     }
     else if (key == '4') {
	std::cout << "offset point cloud" << std::endl;
 	viewer.data().point_size = 10;
	viewer.data().clear();	
	viewer.data().set_points(closest_points, Eigen::RowVector3d(1,1,1));
     }
   /* else if (key == '5') {
	std::cout << "igl offset implementation" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V_igl, F_igl);
	viewer.core.align_camera_center(V_igl, F_igl);
     } */
     
     return false;
}

int main(int argc, char *argv[])
{
  // Load mesh with desired distance and resolution
  // argv[1] mesh
  // argv[2] distance away from original mesh
  // argv[3] grid resolution
  // argv[4] radius of influence

  igl::read_triangle_mesh(argv[1], V, F);
  sigma = atof(argv[2]);
  int res = atol(argv[3]);
  double r = atof(argv[4]);
  std::cout << "distance: " << sigma << "resolution: " << res << std::endl;
  marching_cubes_offset(V, F, sigma, res, r, V_mc, F_mc, V_mcr, F_mcr, closest_points);
  std::cout << "marching cubes: " << std::endl;
  validate(V, F, V_mc, F_mc, sigma, int_dist_o);

  std::cout << "marching cubes rf: " << std::endl;
  validate(V, F, V_mcr, F_mcr, sigma, int_dist_n);

/*  Eigen::MatrixXd GV;
  Eigen::VectorXi side;
  Eigen::VectorXd S;
  igl::copyleft::offset_surface(V, F, sigma, res, igl::SIGNED_DISTANCE_TYPE_DEFAULT, V_igl, F_igl, GV, side, S); */

  // Render the original mesh and the offset mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);
  viewer.launch();
  return EXIT_SUCCESS;
}

