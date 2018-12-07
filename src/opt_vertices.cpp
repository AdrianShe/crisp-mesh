#include "opt_vertices.h"
#include "igl/signed_distance.h"
#include <iostream>

void opt_vertices(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const double sigma, const double lambda,
 Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2,  Eigen::MatrixXd & V, Eigen::MatrixXi & F) {
	
	V.resize(V_2.rows(), 3);
	F.resize(F_2.rows(), 3);
	V_last.resize(V_2.rows(), 3);
	F_last.resize(F_2.rows(), 3);
	V = V_2;
 	F = F_2;
	V_last = V_2;
	F_last = F_2;

	// Do an iteration of gradient descient. 
	Eigen::VectorXd dist;
	Eigen::MatrixXd V_temp;
    	Eigen::VectorXd I;
     	Eigen::MatrixXd closest_point; 
     	Eigen::MatrixXd N;
	Eigen::VectorXd dist_scaled;
	double cur_cost = std::numeric_limits<double>::max();
	double new_cost = 0;

     	igl::signed_distance(V, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, I, closest_point, N);
	Eigen::MatrixXd grad_E;
	grad_E.resize(V.rows(), 3);
	V_temp.resize(V.rows(), 3);
	V_temp = V;
	dist_scaled = dist - sigma * (Eigen::VectorXd::Ones(V.rows()));
	new_cost = std::pow(dist_scaled.norm(), 2);
	std::cout << "Cost: " << new_cost << std::endl;
	while (new_cost < cur_cost) {
		cur_cost = new_cost;	 
		V = V_temp;
 		for (int i = 0; i < V.rows(); i++) {
			Eigen::MatrixXd cur_norm = V.row(i) - closest_point.row(i);
			grad_E.row(i) = 2.0 * dist_scaled(i) * (cur_norm / cur_norm.norm());
		}
		V_temp = V - lambda * grad_E; 
     		igl::signed_distance(V_temp, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), 		std::numeric_limits<double>::max(), dist, I, closest_point, N);
		Eigen::VectorXd dist_scaled = dist - sigma * (Eigen::VectorXd::Ones(V.rows()));
		new_cost = std::pow(dist_scaled.norm(), 2);
		std::cout << "Cost: " << new_cost << std::endl;
	}
}

