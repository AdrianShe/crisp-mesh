#include "validate.h"
#include "igl/signed_distance.h"
#include <iostream>

void validate(const Eigen::MatrixXd & V_1, 
const Eigen::MatrixXi & F_1, 
const Eigen::MatrixXd & V_2, 
const Eigen::MatrixXi & F_2, double sigma) {
	std::cout << "The generated offset mesh has " << V_2.rows() << " vertices and " << F_2.rows() << " faces." << std::endl;
	Eigen::VectorXd dist;
     	Eigen::VectorXd I;
     	Eigen::MatrixXd closest_point; 
     	Eigen::MatrixXd N;
     	igl::signed_distance(V_2, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, I, closest_point, N);
	double min = dist.minCoeff();
	double max = dist.maxCoeff();
	double avg = dist.mean();
	Eigen::ArrayXd sq_dev = dist.array();
        sq_dev = sq_dev - sigma;
	sq_dev = sq_dev * sq_dev;
	double stdev = std::sqrt(sq_dev.mean()); 
	std::cout << "The minimum offset to original distance " << min << std::endl;
	std::cout << "The maximum offset to original distance " << max << std::endl;
	std::cout << "The mean offset to original distance " << avg << std::endl;
	std::cout << "The standard deviation of offset to original distance " << stdev << std::endl;
}
