#ifndef DUAL_CONTOUR
#define DUAL_CONTOUR
#include <Eigen/Core>
// Dual contour implementation on a uniform grid with root finding.
// Based on existing implementation: https://github.com/emilk/Dual-Contouring
// 
// Inputs:
//  dist: isovalues at all grid points
//  grid_pos: grid positions
//  sidex: grid side length x
//  sidey: grid side length y
//  sidez: grid side length z
//  V: vertices of original mesh
//  F: faces of original mesh
//  sigma: offset
//
//  Outputs:
//  (V2, F2) is the offset surface generated from (V, F) 
void dual_contour(const Eigen::VectorXd& dist, const Eigen::MatrixXd& grid_pos,
    const int& sidex, const int& sidey, const int& sidez, 
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double& sigma,
    Eigen::MatrixXd& V2, Eigen::MatrixXi& F2);
#endif
