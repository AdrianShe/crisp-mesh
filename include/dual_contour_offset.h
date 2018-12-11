#ifndef DUAL_CONTOUR_OFFSET
#define DUAL_CONTOUR_OFFSET
#include <Eigen/Core>
// Creates an offset mesh using our dual contouring implementation
// (uniform grid). Calls the implementation function 'dual_contour'.
//  
// Inputs:
//  V_1 mesh vertices
//  F_1 mesh faces
//  sigma offset distance
//  res desired mesh resolution (of the smallest grid side length)
//
//  Outputs:
//  (V_2, F_2) is the offset surface generated from (V_1, F_1) 
void dual_contour_offset(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, 
    const double sigma, const int res, 
    Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2);
#endif
