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
     	igl::signed_distance(V_2, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, I, closest_point, N);
	// double min = dist.minCoeff();
	// double max = dist.maxCoeff();
	double avg = dist.mean();
	Eigen::ArrayXd sq_dev = dist.array();
        sq_dev = sq_dev - sigma;
	sq_dev = sq_dev * sq_dev;
	double stdev = std::sqrt(sq_dev.mean()); 
	// std::cout << "The minimum offset to original distance " << min << std::endl;
	// std::cout << "The maximum offset to original distance " << max << std::endl;
	std::cout << "The mean offset to original distance " << avg << "." << std::endl ;
	std::cout << "The standard deviation of offset to original distance " << stdev << "." << std::endl;

	// Now compute integrated distances over the mesh using quadrature rule on page 422-423 http://www.techmat.vgtu.lt/~inga/Files/Quarteroni-SkaitMetod.pdf
	int_dist.resize(F_2.rows());
	Eigen::VectorXd areas;
	igl::doublearea(V_2, F_2, areas);
	double cum_int = 0.0;
	for (int i = 0; i < F_2.rows(); i++) {
		Eigen::MatrixXd query(7, 3);
		query.row(0) = V_2.row(F_2(i, 0)); // vertices
		query.row(1) = V_2.row(F_2(i, 1)); 
		query.row(2) = V_2.row(F_2(i, 2));
		query.row(3) = (query.row(0) + query.row(1)) / 2; // edge midpoints
 		query.row(4) = (query.row(1) + query.row(2)) / 2; 
		query.row(5) = (query.row(0) + query.row(2)) / 2; 
		query.row(6) = (query.row(3) + query.row(4) + query.row(5)) / 3; // centroid
		dist.setZero();
	     	igl::signed_distance(query, V_1, F_1, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, I, closest_point, N);
		int_dist(i) = (1.0 / 60.0) * (3.0 * (dist(0) + dist(1) + dist(2)) + 8.0 * (dist(3) + dist(4) + dist(5)) + 27.0 * dist(6));
		cum_int = cum_int + areas(i) * int_dist(i); 
		// std::cout << cum_int << std::endl;
	}
	std::cout << "The integrated distance of the offset mesh is " << (cum_int) / (areas.sum()) << "." << std::endl;
Eigen::ArrayXd int_dist_a = int_dist.array() - sigma;
	int_dist = int_dist_a.matrix();
}	
