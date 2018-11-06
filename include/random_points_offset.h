#ifndef RANDOM_PTS_OFFSET
#define RANDOM_PTS_OFFSET
#include <Eigen/Core>
// Sample random points approximately distance sigma from the mesh
//
// Inputs:
//  V_1 mesh vertices
//  F_1 mesh faces
//  sigma offset distance
//  n number of sample points
//
//  Outputs:
//  P #P by 3 array of sample points
//  N_p #Associated surface normals
void random_points_offset(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const double sigma, const int n, Eigen::MatrixXd & P, Eigen::MatrixXd & N_p);
#endif
