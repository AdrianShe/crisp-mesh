#include "validate.h"
#include "igl/signed_distance.h"
#include "igl/doublearea.h"
#include <iostream>

void validate(const Eigen::MatrixXd & V_1, 
const Eigen::MatrixXi & F_1, 
const Eigen::MatrixXd & V_2, 
const Eigen::MatrixXi & F_2, double sigma, Eigen::VectorXd & int_dist) {
	std::cout << "The generated offset mesh has " << V_2.rows() << " vertices and " << F_2.rows() << " faces." << std::endl;
	Eigen::VectorXd dist;
     	Eigen::VectorXd I;
     	Eigen::MatrixXd closest_point; 
     	Eigen::MatrixXd N;

	// Compute average per-vertex error     	
	igl::signed_distance(V_2, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, I, closest_point, N);
	double avg = dist.mean();
	Eigen::ArrayXd sq_dev = dist.array();
        sq_dev = sq_dev - avg;
	sq_dev = sq_dev * sq_dev;
	double stdev = std::sqrt(sq_dev.mean()); 
	std::cout << "Error in mean vertex distance: " << avg - sigma << "." << std::endl ;
	std::cout << "Standard Deviation in Distance: " << stdev << "." << std::endl;
	
	// Compute integrated per-vertex error
	Eigen::MatrixXd query(7 * F_2.rows(), 3);
	for (int i = 0; i <  F_2.rows(); i++) {
		query.row(7 * i) = V_2.row(F_2(i,0));
		query.row(7 * i + 1) = V_2.row(F_2(i,1));
		query.row(7 * i + 2) = V_2.row(F_2(i,2));
		query.row(7 * i + 3) = (V_2.row(F_2(i,0)) + V_2.row(F_2(i,1))) / 2.0;
		query.row(7 * i + 4) = (V_2.row(F_2(i,1)) + V_2.row(F_2(i,2))) / 2.0;
		query.row(7 * i + 5) = (V_2.row(F_2(i,0)) + V_2.row(F_2(i,2))) / 2.0;
		query.row(7 * i + 6) = (V_2.row(F_2(i,0)) + V_2.row(F_2(i,1)) + V_2.row(F_2(i,2))) / 3.0;
	}
	// std::cout << "query points set" << "." << std::endl;

	// Now compute integrated distances over the mesh using quadrature rule on page 422-423 http://www.techmat.vgtu.lt/~inga/Files/Quarteroni-SkaitMetod.pdf
	int_dist.resize(F_2.rows());
	Eigen::VectorXd areas;
	igl::doublearea(V_2, F_2, areas);
	double cum_int = 0.0;
	igl::signed_distance(query, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, I, closest_point, N);
	for (int i = 0; i < F_2.rows(); i++) {
		int_dist(i) = (1.0 / 60.0) * (3.0 * (dist(7 * i) + dist(7 * i + 1) + dist(7 * i + 2)) + 8.0 * (dist(7 * i + 3) + dist(7 * i + 4) + dist(7 * i + 5)) + 27.0 * dist(7 * i + 6));
		cum_int = cum_int + areas(i) * int_dist(i); 
	}
	cum_int = cum_int / areas.sum();
	std::cout << "Error in integrated distance: " << cum_int - sigma << "." << std::endl;
	sq_dev = int_dist.array();
        sq_dev = sq_dev - sq_dev.mean();
	sq_dev = sq_dev * sq_dev;
	stdev = std::sqrt(sq_dev.mean()); 
	std::cout << "Standard Deviation in Distance: " << stdev << "." << std::endl;
	Eigen::ArrayXd int_dist_a = int_dist.array() - sigma;
	int_dist = int_dist_a.matrix();
}	
