#include "opt_vertices_midpoint.h"
#include "igl/signed_distance.h"
#include "igl/doublearea.h"
// #include "igl/copyleft/cgal/remesh_self_intersections.h"
#include <iostream>

void compute_midpoints(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, Eigen::MatrixXd & midpoints) {
	midpoints.resize(3 * F.rows(), 3);	
	for (int i = 0; i < F.rows(); i++) {
		midpoints.row(3 * i) = (V.row(F(i,0)) + V.row(F(i,1))) / 2.0;
		midpoints.row(3 * i + 1) = (V.row(F(i,1)) + V.row(F(i,2))) / 2.0;
		midpoints.row(3 * i + 2) = (V.row(F(i,0)) + V.row(F(i,2))) / 2.0;
	}
}

void compute_normals(const Eigen::MatrixXd & V, const Eigen::MatrixXd & cp, Eigen::MatrixXd & N) {
	N.resize(V.rows(), 3);
	for (int i = 0; i < cp.rows(); i++) {
		N.row(i) = V.row(i) - cp.row(i);
		N.row(i).normalize();
	}	
}

double compute_integral_midpoint( const Eigen::VectorXd & dist, const Eigen::MatrixXi & F, const Eigen::VectorXd & A) {
	double val = 0;
	for (int i = 0; i < F.rows(); i++) {
		val += (A(i) / 3.0) * (std::pow(dist(3 * i), 2) + std::pow(dist(3 * i + 1), 2) + std::pow(dist(3 * i + 2), 2)); 
	}
	return val;
}	


void opt_vertices_midpoint(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const double sigma, const double lambda, const double tol, Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2,  Eigen::MatrixXd & V) {
	V.resize(V_2.rows(), 3);
	// F.resize(F_2.rows(), 3);
	V = V_2;
 	// F = F_2;
	// double tol = 1e-6;
	int num_its = 0;

	// Do an iteration of gradient descent. 

	// Compute current cost
	Eigen::VectorXd dist;
	Eigen::VectorXd dist_scaled;
    	Eigen::VectorXd I;
     	Eigen::MatrixXd closest_point; 
     	Eigen::MatrixXd N;
	double cur_cost = std::numeric_limits<double>::max();
	double new_cost = 0;
	// int min_its_to_do = 1000;

	Eigen::VectorXd A;
	igl::doublearea(V, F_2, A);
	Eigen::MatrixXd midpoints;
	compute_midpoints(V, F_2, midpoints);
	igl::signed_distance(midpoints, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, I, closest_point, N);
	dist_scaled = dist - sigma * (Eigen::VectorXd::Ones(midpoints.rows()));
	new_cost = compute_integral_midpoint(dist_scaled, F_2, A);
	// std::cout << "Initialized" << std::endl;
	std::cout << "Initial Cost: " << new_cost << std::endl;

	// Initialize a bunch of stuff for gradient descent.
	Eigen::MatrixXd grad_E;	
	grad_E.resize(V.rows(), 3);
	Eigen::MatrixXd V_temp;
	V_temp.resize(V.rows(), 3);
	V_temp = V;

	while ( ((new_cost < cur_cost) && std::abs(new_cost - cur_cost) >= tol)) {
		cur_cost = new_cost;	 
		V = V_temp;
 		compute_normals(midpoints, closest_point, N);
		// std::cout << "normal" << std::endl;
		grad_E.setZero(); 
		for (int i = 0; i < F_2.rows(); i++) {
			grad_E.row(F_2(i,0)) += (A(i)/3.0) * (dist_scaled(3 * i) * N.row(3 * i));
			grad_E.row(F_2(i,1)) += (A(i)/3.0) * (dist_scaled(3 * i) * N.row(3 * i));
			grad_E.row(F_2(i,1)) += (A(i)/3.0) * (dist_scaled(3 * i + 1) * N.row(3 * i + 1));
			grad_E.row(F_2(i,2)) += (A(i)/3.0) * (dist_scaled(3 * i + 1) * N.row(3 * i + 1));
			grad_E.row(F_2(i,2)) += (A(i)/3.0) * (dist_scaled(3 * i + 2) * N.row(3 * i + 2));
			grad_E.row(F_2(i,0)) += (A(i)/3.0) * (dist_scaled(3 * i + 2) * N.row(3 * i +  2 )); 
		}
		V_temp = V - lambda * grad_E; 
		// std::cout << "grad" << std::endl;
		compute_midpoints(V_temp, F_2, midpoints);
     		igl::signed_distance(midpoints, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), 		std::numeric_limits<double>::max(), dist, I, closest_point, N);
		dist_scaled = dist - sigma * (Eigen::VectorXd::Ones(midpoints.rows()));
		igl::doublearea(V_temp, F_2, A);
		new_cost = compute_integral_midpoint(dist_scaled, F_2, A);
		// std::cout << "Cost: " << new_cost << " " << std::abs(new_cost - cur_cost) << std::endl;
		num_its++;
		// std::cout << num_its << std::endl;
	}
	std::cout << "Number of iterations: " << num_its << std::endl;
	std::cout << "Final cost: " << cur_cost << std::endl;
	/* Eigen::MatrixXd V_new;
	Eigen::MatrixXi F_new;
	igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
	Eigen::MatrixXd If;
	Eigen::VectorXi J;
	Eigen::VectorXi IM;
	igl::copyleft::cgal::remesh_self_intersections(V,F,params, V_new, F_new, If, J, IM);
	V = V_new;
	F = F_new; */ 
}


