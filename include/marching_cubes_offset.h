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
//  V_2 offset mesh vertices
//  F_2 offset mesh faces 
void marching_cubes_offset(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const double sigma, const int res, Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2);
#endif
