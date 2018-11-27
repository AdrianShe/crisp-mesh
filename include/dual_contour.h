#ifndef DUAL_CONTOUR
#define DUAL_CONTOUR
#include <Eigen/Core>
// Dual contour
// 
// Inputs:
//  TBA
void dual_contour(const Eigen::VectorXd& dist, const Eigen::MatrixXd& grid_pos,
    const int& sidex, const int& sidey, const int& sidez, 
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double& sigma,
    Eigen::MatrixXd& V2, Eigen::MatrixXi& F2);
#endif
