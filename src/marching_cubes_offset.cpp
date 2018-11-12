#include "marching_cubes_offset.h"
#include <igl/copyleft/marching_cubes.h>
#include <igl/copyleft/marching_cubes_root_finding.h>
#include "grid.h"
#include "igl/point_mesh_squared_distance.h"
#include <iostream>

void marching_cubes_offset(
const Eigen::MatrixXd & V_1, 
const Eigen::MatrixXi & F_1, 
const double sigma,  
const int res,
Eigen::MatrixXd & V_2, 
Eigen::MatrixXi & F_2,
Eigen::MatrixXd & V_3,
Eigen::MatrixXi & F_3,
Eigen::MatrixXd & closest_point_cloud)
 {
     // Determine coordinates of bounding box for the offset surface
     double min_x = V_1.col(0).minCoeff() - 2.0 * sigma;
     double min_y = V_1.col(1).minCoeff() - 2.0 * sigma;
     double min_z = V_1.col(2).minCoeff() - 2.0 * sigma;
     double max_x = V_1.col(0).maxCoeff() + 2.0 * sigma;
     double max_y = V_1.col(1).maxCoeff() + 2.0 * sigma;
     double max_z = V_1.col(2).maxCoeff() + 2.0 * sigma;
     Eigen::RowVector3d lower(min_x, min_y, min_z);
     Eigen::RowVector3d upper(max_x, max_y, max_z);
     std::cout << "lower corner coordinates: " << lower << std::endl;
     std::cout << "upper corner coordinates: " << upper << std::endl;
     Eigen::RowVector3d side_lens = upper - lower;
     
     // Scale grid such that min_side_len gets res query points and create it
     double h = side_lens.minCoeff() / res; 
     Eigen::RowVector3d normalized_lens = side_lens / h;
     Eigen::RowVector3i side;
     for (int i = 0; i <= 2; i++) {
	side[i] = ceil(normalized_lens[i]);
     }
     Eigen::MatrixXd grid_pos;
     grid(h, side[0], side[1], side[2], lower, grid_pos);
 
     // Compute grid_pos to mesh distances
     Eigen::VectorXd dist;
     Eigen::VectorXd I;
     Eigen::MatrixXd closest_point; 
     igl::point_mesh_squared_distance(grid_pos, V_1, F_1, dist, I, closest_point);

     // produce mesh using marching cubes from distance data with iso level sigma^2 (as squared distances were returned).
     igl::copyleft::marching_cubes(dist, grid_pos, side[0], side[1], side[2], sigma * sigma, V_2, F_2);
     std::cout << "done normal one" << std::endl;
     // Marching cubes with root finding
     // Sample points close enough to original mesh (V, F)
     int n = 0;
     std::vector<int> close_point_idx;
     double eps = sigma / 5; 
     for (int i = 0; i < grid_pos.rows(); i++) {
	if (std::abs(dist[i] - sigma) <= eps) {
		close_point_idx.push_back(i);
	} 
     }	
     closest_point_cloud.resize(close_point_idx.size(), 3);
  /*   std::uniform_real_distribution<double> x_pos(min_x, max_x);
     std::uniform_real_distribution<double> y_pos(min_y, max_y);
     std::uniform_real_distribution<double> z_pos(min_z, max_z); */
     std::cout << closest_point_cloud.rows() << " " << closest_point_cloud.cols() << std::endl;
     std::cout << close_point_idx.size() << " " << closest_point_cloud.cols() << std::endl;
     for (int i = 0; i < close_point_idx.size(); i++) {
	closest_point_cloud.row(i) = grid_pos.row(close_point_idx[i]);
	// std::cout << i << std::endl;
     }
  //  std::cout << closest_point_cloud << std::endl;
     igl::copyleft::marching_cubes_root_finding(4 * sigma, -1, closest_point_cloud, grid_pos, side[0], side[1], side[2], igl::copyleft::LOCAL_IMPLICIT_FUNCTION_DEFAULT, V_3, F_3);
} 
