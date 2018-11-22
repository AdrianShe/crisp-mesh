#include "marching_cubes_offset.h"
#include <igl/copyleft/marching_cubes.h>
#include <igl/copyleft/marching_cubes_root_finding.h>
#include "grid.h"
#include "igl/point_mesh_squared_distance.h"
#include "igl/signed_distance.h"
#include <iostream>
#include <random>

void marching_cubes_offset(
const Eigen::MatrixXd  & V_1, 
const Eigen::MatrixXi  & F_1, 
const double sigma,  
const int res,
const double r,
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
     Eigen::MatrixXd N;
     igl::signed_distance(grid_pos, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, sigma - 3 * h, sigma + 3 * h, dist, I, closest_point, N);

     // produce mesh using marching cubes from distance data with iso level sigma 
     igl::copyleft::marching_cubes(dist, grid_pos, side[0], side[1], side[2], sigma, V_2, F_2);

     // Marching cubes with root finding
     // Sample points close enough to original mesh (V, F)
  /*   std::default_random_engine gen;
     std::uniform_real_distribution<double> x_pos(min_x, max_x);
     std::uniform_real_distribution<double> y_pos(min_y, max_y);
     std::uniform_real_distribution<double> z_pos(min_z, max_z);
     Eigen::MatrixXd sample_point_cloud(100000, 3);
     for (int i = 0; i < sample_point_cloud.rows(); i++) {
	Eigen::RowVector3d point(x_pos(gen), y_pos(gen), z_pos(gen));
	sample_point_cloud.row(i) = point;
     }	
     igl::signed_distance(sample_point_cloud, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT,  sigma - 3 * h, sigma + 3 * h, dist, I, closest_point, N); 



     // compute number of points in sample_point_cloud lying distance sigma from the mesh
     int n = 0;
     for (int i = 0; i < sample_point_cloud.rows(); i++) {
	if (dist[i] <= sigma) {
		n++;
	}
     }
     closest_point_cloud.resize(n, 3);
     n = 0;
     for (int i = 0; i < sample_point_cloud.rows(); i++) {
	if (dist[i] <= sigma) {
		closest_point_cloud.row(n) = sample_point_cloud.row(i);
		n++;
	}
     } */
     igl::copyleft::marching_cubes_root_finding(V_1, F_1, r, sigma, V_1, grid_pos, (unsigned int) side[0], (unsigned int) side[1], (unsigned int) side[2], igl::copyleft::LocalImplicitFunction::SIGNED_DISTANCE, V_3, F_3);
} 
