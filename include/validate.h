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
void validate(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const Eigen::MatrixXd & V_2, const Eigen::MatrixXi & F_2, double sigma);
#endif
