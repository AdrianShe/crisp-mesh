#include "grid.h"
#include <iostream>

void grid(const double h, 
const int nx, 
const int ny, 
const int nz, 
const Eigen::RowVector3d lower, 
Eigen::MatrixXd & grid_pos) {
     grid_pos.resize(nx * ny * nz, 3); 
     for (int i = 0; i < nx; i++) {
	for (int j = 0; j < ny; j++) {
		for (int k = 0; k < nz; k++) {
			grid_pos.row(i + nx * (j + ny * k)) = lower + h * Eigen::RowVector3d(i, j, k);
		}
	}
     }
}


