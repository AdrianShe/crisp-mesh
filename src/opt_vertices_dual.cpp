#include "opt_vertices_dual.h"
#include "igl/signed_distance.h"
#include "igl/vertex_triangle_adjacency.h"
#include "igl/doublearea.h"
#include <iostream>

double compute_cost(const double& sigma, const Eigen::VectorXd& A, 
    const Eigen::MatrixXd& proj_C, const Eigen::MatrixXd& N, 
    const Eigen::MatrixXd& V, const std::vector<std::vector<int>>& VF) {
  double output = 0;
  for (int i = 0; i < V.rows(); ++i) {
    for (int j = 0; j < VF[i].size(); ++j) {
      int f = VF[i][j];
      double error = (V.row(i) - proj_C.row(f)).dot(N.row(f)) - sigma;
      output += A(f)*error*error;
    } // end loop j
  } // end loop i
  return 0.5*output;
}

void compute_centroids(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, 
    Eigen::MatrixXd& C) {

  C.resize(F.rows(), 3);
  for (int i = 0; i < C.rows(); ++i) {
    C.row(i) = V.row(F(i,0)) + V.row(F(i,1)) + V.row(F(i,2));
    C.row(i) /= 3.0;
  } // end loop i

}

void approximate_normals(const Eigen::MatrixXd& C, const Eigen::MatrixXd& proj_C,
    Eigen::MatrixXd& N) {

  N = C - proj_C;
  for (int i = 0; i < N.rows(); ++i)
    N.row(i).normalize();
}

void opt_vertices_dual(const Eigen::MatrixXd & V_1, const Eigen::MatrixXi & F_1, 
    const double sigma, const double lambda,const double tol,
    Eigen::MatrixXd & V_2, Eigen::MatrixXi & F_2,  Eigen::MatrixXd & V) {
	
  bool area_flag = false;

  // Initialize V
	V = V_2;

  // Other userful variables
	Eigen::MatrixXd V_temp, C;
  V_temp = V;
	Eigen::VectorXd I, A, dist;
  A = Eigen::VectorXd::Ones(F_2.rows());
 	Eigen::MatrixXd closest_point, N; 

	double cur_cost = std::numeric_limits<double>::max();
	double new_cost = 0;

  // Determine neighbouring faces
  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(V_2, F_2, VF, VFi);

  // Initial cost
  compute_centroids(V, F_2, C);
  igl::signed_distance(C, V_1, F_1, 
      igl::SIGNED_DISTANCE_TYPE_DEFAULT, 
      std::numeric_limits<double>::min(), 
      std::numeric_limits<double>::max(), 
      dist, I, closest_point, N);
  igl::doublearea(V, F_2, A);
  A *= 0.5;
  approximate_normals(C, closest_point, N);

  new_cost = compute_cost(sigma, A, closest_point, N, V, VF);
	std::cout << "Initial Cost: " << new_cost << std::endl;

	
	int num_its = 0;
	// int min_its_to_do = 1000; 

  // loop while cost decreases
	while  ((new_cost < cur_cost) && std::abs(new_cost - cur_cost) >= tol)  {
		cur_cost = new_cost;	 
		V = V_temp;

    // Compute update on V_temp
    compute_centroids(V, F_2, C);
    igl::signed_distance(C, V_1, F_1, 
      igl::SIGNED_DISTANCE_TYPE_DEFAULT, 
      std::numeric_limits<double>::min(), 
      std::numeric_limits<double>::max(), 
      dist, I, closest_point, N);
    igl::doublearea(V, F_2, A);
    A *= 0.5;
    approximate_normals(C, closest_point, N);

    for (int i = 0; i < V_temp.rows(); ++i) {
      for (int j = 0; j < VF[i].size(); ++j) {
        int f = VF[i][j];
        V_temp.row(i) -= lambda*A(f)*(
            N.row(f).dot(V.row(i) - closest_point.row(f)) - sigma
          )*N.row(f);
      } // end loop j
    } // end loop i

    // Compute new cost
	  new_cost = compute_cost(sigma, A, closest_point, N, V_temp, VF);
	//std::cout << "Cost: " << new_cost << std::endl;
	num_its++;
	}  // end while loop
	std::cout << "Number of iterations: " << num_its << std::endl;
	std::cout << "Final cost: " << cur_cost << std::endl;
  // std::cout << "Cost: " << new_cost << std::endl;
}
