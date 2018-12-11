#ifndef VALIDATE
#define VALIDATE
#include <Eigen/Dense>
// Test the quality of a purported offset mesh (V_2, F_2) of distance sigma from (V_1, F_1)
// 
// Inputs:
//  V_1, F_1 original mesh
//  V_2, F_2 purported offset mesh
//  sigma purported distance sigma
// 
// Outputs:
//  prints error in average distance of vertices in V_2 to (V_1, F_1) and its standard deviation
//  prints error in integrated distance of (V_2, F_2) 
//  returns int_dist, the integrated signed distance field over each triangle in the mesh (normalized by area) and its standard deviation.
void validate(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const Eigen::MatrixXd & V_2, const Eigen::MatrixXi & F_2, double sigma, Eigen::VectorXd & int_dist);
#endif
