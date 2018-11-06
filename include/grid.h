#ifndef GRID_H
#define GRID_H
#include <Eigen/Core>
// Creates a uniform cubic grid 
// 
// Inputs:
//  h grid spacing
//  nx no. of x points
//  ny no. of y points
//  nz no. of z points
//  lower the bottommost frontmost leftmost point on the grid
// 
//  Outputs:
//  grid the grid as a (nx * ny * nz) x 3 matrix. 
void grid(const double h, const int nx, const int ny, const int nz, const Eigen::RowVector3d lower, Eigen::MatrixXd & grid_pos);
#endif
