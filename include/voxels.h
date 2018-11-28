#ifndef VOXELS
#define VOXELS
#include <Eigen/Core>
#include <iostream>
#include <vector>

//  TBA
namespace crisp {

class Voxels {
 public:
  Voxels();
  Voxels(const int& sidex, const int& sidey, const int& sidez) 
  : side_{sidex, sidey, sidez} {
    voxels_.resize(sidex*sidey*sidez);
    local_edges_ = Eigen::MatrixXi(12, 2);
    local_edges_ <<
        0,1, 0,2, 0,4, 1,3, 1,5, 2,3, 2,6, 3,7, 4,5, 4,6, 5,7, 6,7;
  }

  // getters
  int getVID(const int& i, const int& j, const int& k) const {
    return voxels_[i + side_[0]*(j + side_[1]*k)].vid;
  }

  Eigen::VectorXi getCorners(
      const int& i, const int& j, const int& k) const {
    Eigen::VectorXi output(8);
    output(0) = coord(i+0,j+0,k+0);
    output(1) = coord(i+0,j+0,k+1);
    output(2) = coord(i+0,j+1,k+0);
    output(3) = coord(i+0,j+1,k+1);
    output(4) = coord(i+1,j+0,k+0);
    output(5) = coord(i+1,j+0,k+1);
    output(6) = coord(i+1,j+1,k+0);
    output(7) = coord(i+1,j+1,k+1);
    return output;
  }

  // Eigen::MatrixXi getFarEdges(
  //   const int& i, const int& j, const int& k) const {
  //   Eigen::MatrixXi output(3,2);
  //   // 3,7
  //   output(0,0) = coord(0,1,1);
  //   output(0,1) = coord(1,1,1);

  //   // 5,7
  //   output(1,0) = coord(1,0,1);
  //   output(1,1) = coord(1,1,1);

  //   // 6,7
  //   output(2,0) = coord(1,1,0);
  //   output(2,1) = coord(1,1,1);
  // }

  int getDim(const int& x) const {
    return side_[x] - 1;
  }

  Eigen::MatrixXi getLocalEdges() const {
    return local_edges_;
  }

  Eigen::MatrixXd getVertices() const {
    Eigen::MatrixXd V(num_vertex_, 3);
    for (int i = 0; i < num_vertex_; ++i) {
      V.row(i) = ordered_vertices_[i];
    }
    return V;
  }

  // setters
  void setVertex(const int& i, const int& j, const int& k, 
      const Eigen::RowVector3d& pos) {
    Voxel& voxel = voxels_[i + side_[0]*(j + side_[1]*k)];
    voxel.pos = pos;
    if (voxel.vid == -1) {
      voxel.vid = num_vertex_;
      ++num_vertex_;
      ordered_vertices_.push_back(pos);
    }
    else {
      ordered_vertices_[voxel.vid] = pos;
    }
  }

  // helper functions
  int coord(const int& i, const int& j, const int& k) const {
    return i + side_[0]*(j + side_[1]*k);
  }

  void printDebug() const{
    std::cout << "side: " << side_[0] << ", " << side_[1] << ", " << side_[2] << std::endl;
  }
 private:
  struct Voxel {
    int vid = -1; // voxel id
    Eigen::RowVector3d pos;
  };

  // Note: we use same dimensions as grid vertices (not correct)
  int side_[3] = {0,0,0}; // voxel grid
  int num_vertex_ = 0;
  std::vector<Voxel> voxels_;
  std::vector<Eigen::RowVector3d> ordered_vertices_; 

  Eigen::MatrixXi local_edges_;
};

} // namespace crisp

#endif
