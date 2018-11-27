#ifndef DUAL_CONTOUR_OFFSET
#define DUAL_CONTOUR_OFFSET
#include <Eigen/Core>
// Outputs grid values as csv
// 
// Inputs:
//  TBA
void dual_contour_offset(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, 
    const double sigma, const int res, 
    Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2);
#endif
