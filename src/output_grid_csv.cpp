#include "output_grid_csv.h"
#include <igl/copyleft/marching_cubes.h>
#include <igl/copyleft/marching_cubes_root_finding.h>
#include "grid.h"
#include <igl/point_mesh_squared_distance.h>
#include <igl/signed_distance.h>

#include <iostream>
#include <random>
#include <igl/flood_fill.h>

void output_grid_csv(
    const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, 
    const double sigma, const int res, 
    Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2) {

  // Determine coordinates of bounding box for the offset surface
  double min_x = V.col(0).minCoeff() - 2.0 * sigma;
  double min_y = V.col(1).minCoeff() - 2.0 * sigma;
  double min_z = V.col(2).minCoeff() - 2.0 * sigma;
  double max_x = V.col(0).maxCoeff() + 2.0 * sigma;
  double max_y = V.col(1).maxCoeff() + 2.0 * sigma;
  double max_z = V.col(2).maxCoeff() + 2.0 * sigma;
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
  igl::signed_distance(grid_pos, V, F, igl::SIGNED_DISTANCE_TYPE_DEFAULT, -sigma - 3 * h, sigma + 3 * h, dist, I, closest_point, N);
  igl::flood_fill(side, dist);


  // produce mesh using marching cubes from distance data with iso level sigma 
  igl::copyleft::marching_cubes(dist, grid_pos, side[0], side[1], side[2], sigma, V_2, F_2);

  Eigen::MatrixXd dist_closest = grid_pos - closest_point;

  std::ofstream myfile;
  myfile.open ("implicit.csv");
  myfile << side[0] << "," << side[1] << "," << side[2] << std::endl;
  for (int i = 0; i < dist.rows(); ++i) {
    // Eigen::Vector3d normal(dist_dx(i), dist_dy(i), dist_dz(i));
    Eigen::Vector3d normal = dist_closest.row(i);
    
    if (normal.norm() > 0.001)
      normal.normalize();

    myfile << dist(i) << "," << normal[0] << "," << normal[1] << "," << normal[2] << std::endl;

  }
} 

