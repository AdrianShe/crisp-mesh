#ifndef OPT_VERTICES_MDPT
#define OPT_VERTICES_MDPT
#include <Eigen/Core>
// Creates a offset mesh by optimizing the vertices
// 
// Inputs:
//  V_1 original vertices
//  F_1 original faces
//  V_2 mesh to improve vertices
//  F_2 mesh to improve vertices
//  sigma offset distance
//
//  Outputs:
//  (V, F) is the offset surface generated from (V_2, F_2) using marching cubes method 
void opt_vertices_midpoint(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const double sigma, const double lambda, const double tol, Eigen::MatrixXd & V_2,  Eigen::MatrixXi & F_2,  Eigen::MatrixXd & V);
#endif
