#include "dual_contour.h"
#include "voxels.h"

#include <igl/signed_distance.h>
#include <igl/pinv.h>
#include <igl/copyleft/quadprog.h>
#include <igl/per_face_normals.h>


// Root finding on an edge by bisection
void bisectionSearch(
    const Eigen::RowVectorXd& p0, const Eigen::RowVectorXd& p1,
    const double& val0, const double& val1,
    const double& sigma, const double& tol,
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
    Eigen::RowVectorXd& c, Eigen::RowVectorXd& n 
  ) {
  Eigen::VectorXi I;
  Eigen::MatrixXd N, closest_point;

  // make sure 0 is the side less than sigma
  Eigen::RowVectorXd q0, q1;
  double qval0, qval1;
  if (val0 < sigma) {
    q0 = p0;
    q1 = p1;
    qval0 = val0;
    qval1 = val1;
  }
  else {
    q0 = p1;
    q1 = p0;
    qval0 = val1;
    qval1 = val0;
  }

  Eigen::RowVectorXd d_10 = q1 - q0;
  double curr = d_10.norm();
  double original_length = curr;
  d_10  /= curr;

  curr *= 0.5;
  c = q0 + curr*d_10;
  Eigen::VectorXd cval;
  igl::signed_distance(c, V, F, igl::SIGNED_DISTANCE_TYPE_DEFAULT, -3*sigma, 3*sigma, cval, I, closest_point, N);
  while (std::fabs(cval(0) - sigma) > tol) {
    if (cval(0) < sigma) {
      // q0 = c;
      curr *= 0.5;
      c += curr*d_10;
    }
    else {
      // q1 = c;
      curr *= 0.5;
      c -= curr*d_10;
    }
    igl::signed_distance(c, V, F, igl::SIGNED_DISTANCE_TYPE_DEFAULT, -3*sigma, 3*sigma, cval, I, closest_point, N);
  }

  double dotp = d_10.dot(c - q0);
  if (dotp < 0 || dotp >= original_length)
    throw std::runtime_error("[bisectionSearch] out of bounds.");

  // c = (q1 + q0).normalized();
  // igl::signed_distance(c, V, F, igl::SIGNED_DISTANCE_TYPE_DEFAULT, -3*sigma, 3*sigma, cval, I, closest_point, N);
  n = (c - closest_point).normalized();

  if (std::fabs(n.norm() - 1.0) > tol)
    throw std::runtime_error("[bisectionSearch] normal not unit.");
  // n(0) = I(0);
}

// dual contour implementation on uniform grid
void dual_contour(const Eigen::VectorXd& dist, const Eigen::MatrixXd& grid_pos,
    const int& sidex, const int& sidey, const int& sidez, 
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double& sigma,
    Eigen::MatrixXd& V2, Eigen::MatrixXi& F2) {

  // Setup useful variables
  crisp::Voxels voxels(sidex, sidey, sidez);
  std::vector<Eigen::RowVector3d> normals;
  std::vector<Eigen::RowVector3d> intersects;
  std::vector<double> dots;

  Eigen::RowVectorXd x, n;
  Eigen::VectorXd y;
  std::vector<Eigen::RowVectorXi> faces;
  Eigen::VectorXd ce;
  Eigen::MatrixXd CE;
  CE.resize(3,0);
  ce.resize(0);
  Eigen::MatrixXd CIT(6,3);
  CIT << 1,0,0, 0,1,0, 0,0,1, -1,0,0, 0,-1,0, 0,0,-1;
  double push_factor = 1e-1;
  double bisection_tol = 1e-3;

  Eigen::MatrixXd face_norms;
  igl::per_face_normals(V, F, face_norms);

  // Construct vertices //

  // loop through all voxels
  for (int i = 0; i < voxels.getDim(0); ++i) {
    for (int j = 0; j < voxels.getDim(1); ++j) {
      for (int k = 0; k < voxels.getDim(2); ++k) {

        bool inside[8];
        int num_inside = 0;
        int num_x = 0;

        Eigen::VectorXi corners = voxels.getCorners(i,j,k);
        for (int ii = 0; ii < corners.rows(); ++ii) {
          inside[ii] = dist(corners(ii)) <= sigma;
          if (inside[ii])
            ++num_inside;
        } // end loop ii

        // skip if fully inside or outside
        if (num_inside == 0 || num_inside == 8)
          continue;

        normals.clear();
        dots.clear();
        const Eigen::MatrixXi edges = voxels.getLocalEdges();
        for (int ii = 0; ii < edges.rows(); ++ii) {
          const Eigen::RowVectorXi& e = edges.row(ii);
            
          // skip if no crossing on edge  
          if (inside[e(0)] == inside[e(1)])
            continue;

          // determine surface-edge intersection
          const Eigen::RowVectorXd& p0 = grid_pos.row(corners(e(0)));
          const Eigen::RowVectorXd& p1 = grid_pos.row(corners(e(1))); 

          // root finding
          bisectionSearch(p0, p1, dist(corners(e(0))), dist(corners(e(1))),
              sigma, bisection_tol, V, F, x, n);


          normals.push_back(n);
          intersects.push_back(x);
          dots.push_back(n.dot(x));
        } // end loop ii

        if (normals.size() < 3)
          throw std::runtime_error("[dual_contour] Less than 3 intersections.");

        // add weight to voxel centre
        Eigen::RowVectorXd x1 = grid_pos.row(corners(0));
        Eigen::RowVectorXd x2 = grid_pos.row(corners(7));
        {
          Eigen::VectorXi I;
          Eigen::MatrixXd N, closest_point;
          Eigen::VectorXd cval;
          x = (x1 + x2)/2.0;      

          normals.push_back( push_factor*Eigen::RowVector3d(1,0,0) );
          dots.push_back(x.dot(normals.back()));
          
          normals.push_back( push_factor*Eigen::RowVector3d(0,1,0) );
          dots.push_back(x.dot(normals.back()));

          normals.push_back( push_factor*Eigen::RowVector3d(0,0,1) );
          dots.push_back(x.dot(normals.back()));

        }

        // solve for vertex
        Eigen::MatrixXd A(normals.size(), 3);
        Eigen::VectorXd b(normals.size());
        for (int ii = 0; ii < normals.size(); ++ii) {
          A.row(ii) = normals[ii];
          b(ii) = dots[ii];
        } // end loop ii

        Eigen::VectorXd ci0(6);
        ci0 << -x1(0), -x1(1), -x1(2), x2(0), x2(1), x2(2);

        bool success = igl::copyleft::quadprog(
            A.transpose()*A, 
            -b.transpose()*A,  
            CE, 
            ce,  
            CIT.transpose(), 
            ci0, 
            y
          );

        if (std::isnan(y(0)) || std::isnan(y(1)) || std::isnan(y(2)) )
          throw std::runtime_error("[dual_contour] nan.");

          voxels.setVertex(i, j, k, y.transpose());
      } // end loop k
    } // end loop j
  } // end loop i

  // Construct faces //

  Eigen::MatrixXi far_edges(3,2);
  far_edges << 3,7, 5,7, 6,7;
  // std::cout << far_edges << std::endl;

  // loop through all voxels
  for (int i = 0; i < voxels.getDim(0); ++i) {
    for (int j = 0; j < voxels.getDim(1); ++j) {
      for (int k = 0; k < voxels.getDim(2); ++k) {

        int v0 = voxels.getVID(i, j, k);
        if (v0 == -1) {
          continue;
        }

        bool inside[8];
        // int num_inside = 0;

        Eigen::VectorXi corners = voxels.getCorners(i,j,k);
        for (int ii = 0; ii < corners.rows(); ++ii) {
          inside[ii] = dist(corners(ii)) <= sigma;
        } // end loop ii


        for (int ii = 0; ii < 3; ++ii) {

          if (inside[far_edges(ii, 0)] == inside[far_edges(ii, 1)])
            continue;

          int v1, v2, v3;

          if (ii == 0) {
            v1 = voxels.getVID(i, j, k+1);
            v2 = voxels.getVID(i, j+1, k);
            v3 = voxels.getVID(i, j+1, k+1);
          }
          else if (ii == 1) {
            v1 = voxels.getVID(i, j, k+1);
            v2 = voxels.getVID(i+1, j, k);
            v3 = voxels.getVID(i+1, j, k+1);
          }
          else {
            v1 = voxels.getVID(i, j+1, k);
            v2 = voxels.getVID(i+1, j, k);
            v3 = voxels.getVID(i+1, j+1, k);
          }

          if (v1 < 0 || v2 < 0 || v3 < 0) {
            throw std::runtime_error("[dual_contour] -1 vertex id.");
            // continue;
          }

          Eigen::RowVector3i t0, t1;
          if (inside[far_edges(ii,0)] != (ii == 1)) {
            t0 << v0, v3, v1;
            t1 << v0, v2, v3;
          }
          else {
            t0 << v0, v1, v3;
            t1 << v0, v3, v2; 
          }
          faces.push_back(t0);
          faces.push_back(t1);

        } // end loop ii

      } // end loop k
    } // end loop j
  } // end loop i

  // get vertices
  V2 = voxels.getVertices();

  // get faces
  F2.resize(faces.size(), 3);
  std::cout << faces.size() << std::endl;
  for (int i = 0; i < faces.size(); ++i) {
    F2.row(i) = faces[i];
  }
} 

