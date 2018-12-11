#ifndef OPT_VERTICES
#define OPT_VERTICES
#include <Eigen/Core>
// Creates a offset mesh by optimizing the vertices using vertex optimization
// 
// Inputs:
//  V_1 original vertices
//  F_1 original faces
//  sigma offset distance
//  lambda gradient descent parameter
//  tol tolerance when to step gradient descent
//  V_2 initial guess for offset surface
//  F_2 initial guess for offset surface
//
//  Outputs:
//  V updated vertices from (V_2, F_2) to minimize the vertex optimization objective function 
void opt_vertices(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const double sigma, const double lambda, const double tol,
  Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2, Eigen::MatrixXd & V);
#endif
