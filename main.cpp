#include "marching_cubes_offset.h"
#include "marching_cubes_offset_rf.h"
#include "opt_vertices.h"
#include "output_grid_csv.h"
#include "dual_contour_offset.h"
#include <igl/write_triangle_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/offset_surface.h>
#include <igl/opengl/glfw/Viewer.h>
#include "igl/jet.h"
#include <Eigen/Core>
#include <string>
#include <iostream>
#include "validate.h"


Eigen::MatrixXd V, V_mc, V_mcr, V_mco, V_igl, V_d, closest_points, C;
Eigen::MatrixXi F, F_mc, F_mcr, F_mco, F_igl, F_d;
Eigen::VectorXd int_dist_o, int_dist_n, int_dist_d;
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
	igl::jet(int_dist_o, - 0.1 * sigma, 0.1 * sigma, C);
	viewer.data().set_colors(C);
     }
     else if (key == '3') {
	std::cout << "offset using marching cubes and optimization" << std::endl;
	viewer.data().clear();
 	viewer.data().set_mesh(V_mco, F_mco);
	igl::jet(int_dist_n, - 0.1 * sigma, 0.1 * sigma, C);
	viewer.data().set_colors(C);
     }
     else if (key == '4') {
	std::cout << "offset using dual contouring" << std::endl;
	viewer.data().clear();	
	viewer.data().set_mesh(V_d, F_d);
	igl::jet(int_dist_d, - 0.1 * sigma, 0.1 * sigma, C);
	viewer.data().set_colors(C);
     }
     return false;
}

int main(int argc, char *argv[])
{
  // Load mesh with desired distance and resolution
  // argv[1] mesh
  // argv[2] distance away from original mesh
  // argv[3] grid resolution


  igl::read_triangle_mesh(argv[1], V, F);
  sigma = atof(argv[2]);
  int res = atol(argv[3]);
  double lambda = atof(argv[4]);

  std::cout << "distance: " << sigma << "resolution: " << res << std::endl;
  marching_cubes_offset(V, F, sigma, res, V_mc, F_mc);
  std::cout << "marching cubes: " << std::endl;
  validate(V, F, V_mc, F_mc, sigma, int_dist_o);

  std::cout << "marching cubes optimization: " <<std::endl;
  opt_vertices(V, F, sigma, lambda, V_mc, F_mc, V_mco, F_mco); 
  validate(V, F, V_mco, F_mco, sigma, int_dist_n);

  // std::cout << "marching cubes rf: " << std::endl;
  // marching_cubes_offset_rf(V, F, sigma, res, V_mcr, F_mcr);
  // validate(V, F, V_mcr, F_mcr, sigma, int_dist_n); 
  // igl::jet(int_dist_n, - 0.1 * sigma, 0.1 * sigma, C);
  // igl::write_triangle_mesh("marching_cubes.off", V_mcr, F_mcr);

  // std::cout << "dual contouring: " <<std::endl;
  // dual_contour_offset(V, F, sigma, res, V_d, F_d);
  // validate(V, F, V_d, F_d, sigma, int_dist_d); 
  // igl::jet(int_dist_d, - 0.1 * sigma, 0.1 * sigma, C);
  // igl::write_triangle_mesh("dual_contouring.off", V_d, F_d);

  // Render the original mesh and the offset mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);
  viewer.launch();
  return EXIT_SUCCESS;
}

