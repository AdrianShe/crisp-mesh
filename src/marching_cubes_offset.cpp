#include "marching_cubes_offset.h"
#include "igl/point_mesh_squared_distance.h"
#include <igl/copyleft/marching_cubes.h>
#include <iostream>

void marching_cubes_offset(
const Eigen::MatrixXd & V_1, 
const Eigen::MatrixXi & F_1, 
const double sigma,  
const int res,
Eigen::MatrixXd & V_2, 
Eigen::MatrixXi & F_2)
 {
     // Create a uniform grid for query points.  
     double min_x = V_1.row(0).minCoeff() - sigma;
     double min_y = V_1.row(1).minCoeff() - sigma;
     double min_z = V_1.row(2).minCoeff() - sigma;
     double max_x = V_1.row(0).maxCoeff() + sigma;
     double max_y = V_1.row(1).maxCoeff() + sigma;
     double max_z = V_1.row(2).maxCoeff() + sigma;
     Eigen::RowVector3d lower(min_x, min_y, min_z);
     Eigen::RowVector3d side_lens(max_x - min_x, max_y - min_y, max_z - min_z);
     
     // Scale grid such that min_side_len gets n query points;
     double h = side_lens.minCoeff() / res; 
     Eigen::RowVector3d normalized_lens = side_lens / h;
     Eigen::Vector3i side;
     for (int i = 0; i <= 2; i++) {
	side[i] = ceil(normalized_lens[i]);
     }

     int nx = side[0];
     int ny = side[1];
     int nz = side[2];
     Eigen::MatrixXd grid_pos(nx * ny * nz, 3); 
     for (int i = 0; i < nx; i++) {
	for (int j = 0; j < ny; j++) {
		for (int k = 0; k < nz; k++) {
			grid_pos.row(i + nx * (j + ny * k)) = lower + Eigen::RowVector3d(h * i, h * j, h * k);
		}
	}
     }
     
     // Compute grid_pos to mesh distances
     Eigen::VectorXd dist;
     Eigen::VectorXd I;
     Eigen::MatrixXd closest_point; 
     igl::point_mesh_squared_distance(grid_pos, V_1, F_1, dist, I, closest_point);

     // produce mesh from marching cubes
     igl::copyleft::marching_cubes(dist, grid_pos, side[0], side[1], side[2], sigma * sigma, V_2, F_2);
}
