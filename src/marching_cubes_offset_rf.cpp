#include "marching_cubes_offset.h"
#include <igl/copyleft/marching_cubes_root_finding.h>
#include "grid.h"
#include "igl/point_mesh_squared_distance.h"
#include "igl/signed_distance.h"
#include <iostream>

void marching_cubes_offset_rf(
const Eigen::MatrixXd  & V_1, 
const Eigen::MatrixXi  & F_1, 
const double sigma,  
const int res,
Eigen::MatrixXd & V_2, 
Eigen::MatrixXi & F_2) {
  // Determine coordinates of bounding box for the offset surface
  double min_x = V_1.col(0).minCoeff() - 2.0 * sigma;
  double min_y = V_1.col(1).minCoeff() - 2.0 * sigma;
  double min_z = V_1.col(2).minCoeff() - 2.0 * sigma;
  double max_x = V_1.col(0).maxCoeff() + 2.0 * sigma;
  double max_y = V_1.col(1).maxCoeff() + 2.0 * sigma;
  double max_z = V_1.col(2).maxCoeff() + 2.0 * sigma;
  Eigen::RowVector3d lower(min_x, min_y, min_z);
  Eigen::RowVector3d upper(max_x, max_y, max_z);
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
  igl::signed_distance(grid_pos, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, 
      std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), 
      dist, I, closest_point, N);

  // produce mesh using marching cubes from distance data with iso level sigma 
  igl::copyleft::marching_cubes_root_finding(V_1, F_1, 0, sigma, V_1, grid_pos, 
      (unsigned int) side[0], (unsigned int) side[1], (unsigned int) side[2], 
      igl::copyleft::LocalImplicitFunction::SIGNED_DISTANCE, V_2, F_2);
} 
