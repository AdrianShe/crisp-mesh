#include "random_points_offset.h"

#include "igl/random_points_on_mesh.h"
#include "igl/per_face_normals.h"
#include "igl/point_mesh_squared_distance.h"

void random_points_offset(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, const double sigma, const int n, Eigen::MatrixXd & P, Eigen::MatrixXd & N_p) {
	Eigen::MatrixXd O_p(n, 3);
        P.resize(2 * n, 3);
 	Eigen::VectorXi FI(n);
	igl::random_points_on_mesh(n, V_1, F_1, O_p, FI);

        // Convert barycentric coordinates into Euclidean coordinates
	for (int i = 0; i < n; i++) {
		Eigen::RowVector3i cur_face = F_1.row(FI(i)); 
		O_p.row(i) = O_p(i, 0) * V_1.row(cur_face(0)) + O_p(i, 1) * V_1.row(cur_face(1)) + O_p(i,2) * V_1.row(cur_face(2)); 
	}

	// For each point offset by distance sigma away from the unit surface normal
	Eigen::MatrixXd N;
	igl::per_face_normals(V_1, F_1, Eigen::RowVector3d::Zero(), N);
	N_p.resize(2 * n, 3);
	for (int i = 0; i < n; i++) {
		P.row(2 * i) = O_p.row(i) + sigma * N.row(FI(i)); 
		P.row(2 * i + 1) = O_p.row(i) - sigma * N.row(FI(i)); 
		N_p.row(2 * i) = N.row(FI(i));
		N_p.row(2 * i + 1) = - N.row(FI(i)); 
	}

	// Compute the average distance from the points to the unit sigma
        Eigen::VectorXd dist;
        Eigen::VectorXd I;
        Eigen::MatrixXd closest_point; 
        igl::point_mesh_squared_distance(P, V_1, F_1, dist, I, closest_point);
        std::cout << "average squared sample-mesh distance: " << dist.sum() / n << std::endl;
}
