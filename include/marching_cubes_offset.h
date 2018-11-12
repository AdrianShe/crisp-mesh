#ifndef MARCHING_CUBES_OFFSET
#define MARCHING_CUBES_OFFSET
#include <Eigen/Core>
// Creates a offset mesh using the marching cubes method
// 
// Inputs:
//  V_1 mesh vertices
//  F_1 mesh faces
//  sigma offset distance
// 
//  Outputs:
//  grid: (nx * ny * nz) x 3 grid of squared distances from (V, F)
void marching_cubes_offset(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const double sigma, const int res, Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2, Eigen::MatrixXd & V_3, Eigen::MatrixXi & F_3, Eigen::MatrixXd & closest_point_cloud);
#endif
