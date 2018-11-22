#ifndef OUTPUT_GRID_CSV
#define OUTPUT_GRID_CSV
#include <Eigen/Core>
// Outputs grid values as csv
// 
// Inputs:
//  TBA
void output_grid_csv(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, 
    const double sigma, const int res, 
    Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2);
#endif
