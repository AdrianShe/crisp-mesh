#include "marching_cubes_offset.h"
#include "marching_cubes_offset_rf.h"
#include "opt_vertices_midpoint.h"
#include "opt_vertices.h"
#include "output_grid_csv.h"
#include "dual_contour_offset.h"
#include "opt_vertices_dual.h"
#include <igl/write_triangle_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/offset_surface.h>
#include <igl/opengl/glfw/Viewer.h>
#include "igl/jet.h"
#include <Eigen/Core>
#include <string>
#include <iostream>
#include "validate.h"


Eigen::MatrixXd V, V_mc, V_mcrf, V_mcvo, V_mcdo, V_mcmo, V_d, C;
Eigen::MatrixXi F, F_mc, F_mcrf, F_d;
Eigen::VectorXd int_dist_o, int_dist_rf, int_dist_vo, int_dist_mo, int_dist_do, int_dist_d;
double sigma;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
     switch(key) {
	case '1':
	 std::cout << "original" << std::endl;
	 viewer.data().clear();
 	 viewer.data().set_mesh(V, F);
	 viewer.core.align_camera_center(V, F);
	 C = Eigen::RowVector3d(0,0,0);
	 break;
	case '2':
	 std::cout << "offset using marching cubes" << std::endl;
	 viewer.data().clear();
 	 viewer.data().set_mesh(V_mc, F_mc);
	 igl::jet(int_dist_o, - 0.1 * sigma, 0.1 * sigma, C);
	 break;
	case '3':
	 std::cout << "offset using marching cubes and vertex optimization" << std::endl;
	 viewer.data().clear();
 	 viewer.data().set_mesh(V_mcvo, F_mc);
	 igl::jet(int_dist_vo, - 0.1 * sigma, 0.1 * sigma, C);
	 break;
	case '4':
	 std::cout << "offset using marching cubes and dual optimization" << std::endl;
	 viewer.data().clear();
 	 viewer.data().set_mesh(V_mcdo, F_mc);
	 igl::jet(int_dist_do, - 0.1 * sigma, 0.1 * sigma, C);
	 break;
	case '5':
	 std::cout << "offset using marching cubes and midpoint optimization" << std::endl;
	 viewer.data().clear();
 	 viewer.data().set_mesh(V_mcmo, F_mc);
	 igl::jet(int_dist_mo, - 0.1 * sigma, 0.1 * sigma, C);
	 break;
	case '6':
 	 std::cout << "offset using marching cubes and root finding" << std::endl;
	 viewer.data().clear();
 	 viewer.data().set_mesh(V_mcrf, F_mcrf);
	 igl::jet(int_dist_rf, - 0.1 * sigma, 0.1 * sigma, C);
	 break;
	case '7':
	 std::cout << "offset using dual contouring" << std::endl;
	 viewer.data().clear();	
	 viewer.data().set_mesh(V_d, F_d);
	 igl::jet(int_dist_d, - 0.1 * sigma, 0.1 * sigma, C);
	 break;
	default:
	  break;
     }
     
     viewer.data().set_colors(C);
     return false;
}

int main(int argc, char *argv[])
{
  // Load mesh with desired distance and resolution
  // argv[1] mesh
  // argv[2] distance away from original mesh
  // argv[3] grid resolution
  // argv[4] lambda for gradient descent: vertex optimization
  // argv[5] lambda for gradient descent: dual optimization
  // argv[6] lambda for gradient descent: midpoint optimization
  // argv[7] tolerance for vertex opt
  // argv[8] tolerance for dual opt
  // argv[9] tolernace for midpoint opt

  igl::read_triangle_mesh(argv[1], V, F);
  sigma = atof(argv[2]);
  int res = atol(argv[3]);
  double lambda_1 = atof(argv[4]);
  double lambda_2 = atof(argv[5]);
  double lambda_3 = atof(argv[6]);
  double tol_1 = atof(argv[7]);
  double tol_2 = atof(argv[8]);
  double tol_3 = atof(argv[9]);

  std::cout << "distance: " << sigma << " resolution: " << res << std::endl;

  std::cout << "marching cubes: " << std::endl;
  marching_cubes_offset(V, F, sigma, res, V_mc, F_mc);
  validate(V, F, V_mc, F_mc, sigma, int_dist_o);
  igl::write_triangle_mesh("marching_cubes.off", V_mc, F_mc);
  std::cout << "\n" <<std::endl;

  std::cout << "marching cubes vertex optimization: " <<std::endl;
  opt_vertices(V, F, sigma, lambda_1, tol_1, V_mc, F_mc, V_mcvo); 
  validate(V, F, V_mcvo, F_mc, sigma, int_dist_vo);
  igl::write_triangle_mesh("marching_cubes_vertex_opt.off", V_mcvo, F_mc);
    std::cout << "\n" <<std::endl;

  std::cout << "marching cubes dual optimization: " <<std::endl;
  opt_vertices_dual(V, F, sigma, lambda_2, tol_2, V_mc, F_mc, V_mcdo); 
  validate(V, F, V_mcdo, F_mc, sigma, int_dist_do);
  igl::write_triangle_mesh("marching_cubes_dual_opt.off", V_mcdo, F_mc);
  std::cout << "\n" <<std::endl;

  std::cout << "marching cubes midpoint optimization: " <<std::endl;
  opt_vertices_midpoint(V, F, sigma, lambda_3, tol_3, V_mc, F_mc, V_mcmo); 
  validate(V, F, V_mcmo, F_mc, sigma, int_dist_mo);
  igl::write_triangle_mesh("marching_cubes_mid_opt.off", V_mcmo, F_mc);
  std::cout << "\n" <<std::endl; 

  std::cout << "marching cubes root finding: " <<std::endl;
  marching_cubes_offset_rf(V, F, sigma, res, V_mcrf, F_mcrf);
  validate(V, F, V_mcrf, F_mcrf, sigma, int_dist_rf);
  igl::write_triangle_mesh("marching_cubes_rf.off", V_mcrf, F_mcrf);
  std::cout << "\n" <<std::endl;

  std::cout << "dual contouring: " <<std::endl;
  dual_contour_offset(V, F, sigma, res, V_d, F_d);
  validate(V, F, V_d, F_d, sigma, int_dist_d); 
  igl::write_triangle_mesh("dual_contouring.off", V_d, F_d);
  std::cout << "\n" <<std::endl; 

  // Render the original mesh and the offset mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);
  viewer.launch();
  return EXIT_SUCCESS;
}

