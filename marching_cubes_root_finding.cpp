/*===========================================================================*\
 *                                                                           *
 *                                IsoEx                                      *
 *        Copyright (C) 2002 by Computer Graphics Group, RWTH Aachen         *
 *                         www.rwth-graphics.de                              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
 \*===========================================================================*/

#include "marching_cubes_root_finding.h"
#include "marching_cubes_tables.h"
#include "../signed_distance.h"
#include <cmath>
#include <unordered_map>
#include <iostream>

extern const int edgeTable[256];
extern const int triTable[256][2][17];
extern const int polyTable[8][16];

struct EdgeKey1
{
  EdgeKey1(unsigned i0, unsigned i1) : i0_(i0), i1_(i1) {}

  bool operator==(const EdgeKey1& _rhs) const
  {
    return i0_ == _rhs.i0_ && i1_ == _rhs.i1_;
  }

  unsigned i0_, i1_;
};

struct EdgeHash1
{
    std::size_t operator()(const EdgeKey1& key) const {
        std::size_t seed = 0;
        seed ^= key.i0_ + 0x9e3779b9 + (seed<<6) + (seed>>2); // Copied from boost::hash_combine
        seed ^= key.i1_ + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return std::hash<std::size_t>()(seed);
    }
};


class MarchingCubesRF
{
  typedef std::unordered_map<EdgeKey1, unsigned, EdgeHash1> MyMap;
  typedef typename MyMap::const_iterator                  MyMapIterator;

public:
static double implicit_local_function_wyvill_1986(double x1, double y1, double z1, double x2, double y2, double z2, double R) {
  double d_sq = pow((x1 - x2), 2.0) + pow((y1 - y2), 2.0) + pow((z1 - z2), 2.0);
  double dbyr_sq = d_sq/pow(R, 2.0);
  double a = -4.0/9.0;
  double b = 17.0/9;
  double c = -22.0/9;
  double local_field = a * pow(dbyr_sq, 3.0) + b * pow(dbyr_sq, 2.0) + c * dbyr_sq + 1.0;
  return local_field;
}

static double implicit_function_arbitrary_point(const Eigen::MatrixXd &  V, const Eigen::MatrixXi & F, const Eigen::MatrixXd &P, const igl::copyleft::LocalImplicitFunction LOCAL_IMPLICIT_FUNCTION_TYPE, const double x, const double y, const double z, const double R) {

  double field_value = 0;
  int influence_points = 0;

  // get all the points in P that are less than R distance to the test point provided here by (x, y, z)
  if ((LOCAL_IMPLICIT_FUNCTION_TYPE == igl::copyleft::LOCAL_IMPLICIT_FUNCTION_DEFAULT) || (LOCAL_IMPLICIT_FUNCTION_TYPE == igl::copyleft::LOCAL_IMPLICIT_FUNCTION_WYVILL_1986)) {
  for (int i = 0; i < P.rows(); i++) {
    double d_sq = pow((x - P(i, 0)), 2.0) + pow((y - P(i, 1)), 2.0) + pow((z - P(i, 2)), 2.0);
    if (d_sq >= pow(R, 2)) {
      continue;
    }
      field_value += implicit_local_function_wyvill_1986(x, y, z, P(i, 0), P(i, 1), P(i, 2), R);
      influence_points += 1;
     // add more else statements if adding more functions
  }
}
  else if (LOCAL_IMPLICIT_FUNCTION_TYPE == igl::copyleft::SIGNED_DISTANCE) {
     Eigen::VectorXd dist;
     Eigen::VectorXd I;
     Eigen::MatrixXd closest_point; 
     Eigen::MatrixXd N;
     Eigen::MatrixXd query_P(1, 3);
     query_P << x, y ,z;
	igl::signed_distance(query_P, V, F, igl::SIGNED_DISTANCE_TYPE_DEFAULT, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), dist, I, closest_point, N);
	return dist[0];
  }

  // if (influence_points != 0)
  //   // normalize the field_value based on the point contributions
  //   field_value = field_value/double(influence_points);

  return field_value;
}

static void implicit_function_cutoff(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const igl::copyleft::LocalImplicitFunction LOCAL_IMPLICIT_FUNCTION_TYPE, const double &R, const Eigen::MatrixXd &P, double & ref) {
  if ((LOCAL_IMPLICIT_FUNCTION_TYPE == igl::copyleft::LOCAL_IMPLICIT_FUNCTION_DEFAULT) || (LOCAL_IMPLICIT_FUNCTION_TYPE == igl::copyleft::LOCAL_IMPLICIT_FUNCTION_WYVILL_1986)) {
    // refer paper -- we need to choose cutoff such that when two bubbles mix, the resulting surface occupies twice the volume
    // calculations not in the paper but done on board 

    double field_val_on_surface = 0.0;
    for (int i =0; i<P.rows(); i++) {
      field_val_on_surface += implicit_function_arbitrary_point(V, F, P, LOCAL_IMPLICIT_FUNCTION_TYPE, P(i,0), P(i,1),P(i,2),R);
    }

    field_val_on_surface = field_val_on_surface/P.rows();


    std::cout<<"The surface field cutoff computed is "<<field_val_on_surface<<std::endl;
    if (field_val_on_surface==1.0)
      std::cout<<"WARNING: FOR THE RADIUS OF INFLUENCE CHOSEN, THE AVG FUNCTION VALUE AT INPUT POINT CLOUD IS 1 (that is, no other point in the neighborhood!). Consider increasing R!"<<std::endl;
    ref = field_val_on_surface;
  } // add more else statements if adding more functions
}

public:
  MarchingCubesRF(const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
		const double R,
                const double cutoff,
                const Eigen::MatrixXd &P,
                const Eigen::MatrixXd &points,
                const unsigned x_res,
                const unsigned y_res,
                const unsigned z_res,
                const igl::copyleft::LocalImplicitFunction LOCAL_IMPLICIT_FUNCTION_TYPE,
                Eigen::MatrixXd &vertices,
                Eigen::MatrixXi &faces)
  {
    assert(P.cols() == 3);
    assert(points.cols() == 3);

    double ref = 0.5; // some initialization. actual val is set below

    if (cutoff < 0)
      implicit_function_cutoff(V, F, LOCAL_IMPLICIT_FUNCTION_TYPE, R, P, ref);
    else
      ref = cutoff;
    

    if(x_res <2 || y_res<2 ||z_res<2)
      return;
    faces.resize(10000,3);
    int num_faces = 0;

    vertices.resize(10000,3);
    int num_vertices = 0;


    unsigned n_cubes  = (x_res-1) * (y_res-1) * (z_res-1);
    assert(unsigned(points.rows()) == x_res * y_res * z_res);

    unsigned int         offsets_[8];
    offsets_[0] = 0;
    offsets_[1] = 1;
    offsets_[2] = 1 + x_res;
    offsets_[3] =     x_res;
    offsets_[4] =             x_res*y_res;
    offsets_[5] = 1         + x_res*y_res;
    offsets_[6] = 1 + x_res + x_res*y_res;
    offsets_[7] =     x_res + x_res*y_res;
    
    double progress = -0.5;

    for (unsigned cube_it =0 ; cube_it < n_cubes; ++cube_it)
    {
      
      // print progress in steps of 0.5 %
      if (cube_it % int(0.005 * n_cubes) == 0) {
        progress += 0.5;
        std::cout<<"Progress: "<< progress << "% " << std::endl;
      }

      unsigned         corner[8];
      int samples[12];
      unsigned char    cubetype(0);
      unsigned int     i;


      // get point indices of corner vertices
      for (i=0; i<8; ++i)
      {
        // get cube coordinates
        unsigned int _idx = cube_it;
        unsigned int X(x_res-1), Y(y_res-1);
        unsigned int x = _idx % X;  _idx /= X;
        unsigned int y = _idx % Y;  _idx /= Y;
        unsigned int z = _idx;

        // transform to point coordinates
        _idx = x + y*x_res + z*x_res*y_res;

        // add offset
        corner[i] = _idx + offsets_[i];
      }

      //std::cout<<"Trivial Cube: "<<cube_it<<std::endl; 
      // std::cout<<"Realtime Progress: "<<double(cube_it/double(n_cubes))*100<<"%"<<std::endl; 

      // determine cube type
      for (i=0; i<8; ++i) {
        double field = implicit_function_arbitrary_point(V, F, P, LOCAL_IMPLICIT_FUNCTION_TYPE, points(corner[i], 0), points(corner[i], 1), points(corner[i], 2), R);
        if (field > ref)
          cubetype |= (1<<i);
      }


      // trivial reject ?
      if (cubetype == 0 || cubetype == 255)
        continue;

      //std::cout<<"Nontrivial Cube: "<<cube_it<<std::endl;

      // compute samples on cube's edges
      if (edgeTable[cubetype]&1)
        samples[0]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[0], corner[1], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&2)
        samples[1]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[1], corner[2], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&4)
        samples[2]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[3], corner[2], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&8)
        samples[3]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[0], corner[3], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&16)
        samples[4]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[4], corner[5], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&32)
        samples[5]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[5], corner[6], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&64)
        samples[6]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[7], corner[6], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&128)
        samples[7]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[4], corner[7], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&256)
        samples[8]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[0], corner[4], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&512)
        samples[9]  = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[1], corner[5], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&1024)
        samples[10] = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[2], corner[6], vertices, num_vertices, edge2vertex);
      if (edgeTable[cubetype]&2048)
        samples[11] = add_vertex(V, F, R, ref, LOCAL_IMPLICIT_FUNCTION_TYPE, P, points, corner[3], corner[7], vertices, num_vertices, edge2vertex);



      // connect samples by triangles
      for (i=0; triTable[cubetype][0][i] != -1; i+=3 )
      {
        num_faces++;
        if (num_faces > faces.rows())
          faces.conservativeResize(faces.rows()+10000, Eigen::NoChange);

        faces.row(num_faces-1) <<
        samples[triTable[cubetype][0][i  ]],
        samples[triTable[cubetype][0][i+1]],
        samples[triTable[cubetype][0][i+2]];

      }

    }

    vertices.conservativeResize(num_vertices, Eigen::NoChange);
    faces.conservativeResize(num_faces, Eigen::NoChange);

  };

  static int  add_vertex(const Eigen::MatrixXd & V, 
					       const Eigen::MatrixXi & F,
					       const double R, 
                                               const double ref,
                                               const igl::copyleft::LocalImplicitFunction LOCAL_IMPLICIT_FUNCTION_TYPE,
                                               const Eigen::MatrixXd &P,
                                               const Eigen::MatrixXd &points,
                                               unsigned int i0,
                                               unsigned int i1,
                                               Eigen::MatrixXd &vertices,
                                               int &num_vertices,
                                               MyMap &edge2vertex)
  {
    // find vertex if it has been computed already
    MyMapIterator it = edge2vertex.find(EdgeKey1(i0, i1));
    if (it != edge2vertex.end())
      return it->second;
    ;

    double tol = 1e-4; // tolerance while root finding. when to stop.

    // generate new vertex
    const Eigen::Matrix<double, 1, 3> & p0 = points.row(i0);
    const Eigen::Matrix<double, 1, 3> & p1 = points.row(i1);

    double s0 = implicit_function_arbitrary_point(V, F, P, LOCAL_IMPLICIT_FUNCTION_TYPE, points(i0, 0), points(i0, 1), points(i0, 2), R);
    double s1 = implicit_function_arbitrary_point(V, F, P, LOCAL_IMPLICIT_FUNCTION_TYPE, points(i1, 0), points(i1, 1), points(i1, 2), R);

    // do root finding using binary search
    Eigen::Matrix<double, 1, 3> p0_bs = points.row(i0);
    Eigen::Matrix<double, 1, 3> p1_bs = points.row(i1);

    Eigen::Matrix<double, 1, 3> midpt;

    // binary search until tol
    while ((p0_bs-p1_bs).squaredNorm() >= tol) { // squaredNorm is faster than norm()
      // find the middle point
      midpt = (p0_bs + p1_bs)/2.0;

      // compute value here
      double val_here = implicit_function_arbitrary_point(V, F, P, LOCAL_IMPLICIT_FUNCTION_TYPE, midpt(0), midpt(1), midpt(2), R);
      double val_p0_bs = implicit_function_arbitrary_point(V, F, P, LOCAL_IMPLICIT_FUNCTION_TYPE, p0_bs(0), p0_bs(1), p0_bs(2), R);

      // check if mid pt is the root
      if (val_here == ref)
        break; // midpt is the answer
      else if ((val_here - ref) * (val_p0_bs - ref) < 0.0)
        p1_bs = midpt;
      else
        p0_bs = midpt;

    }
    // midpt holds the root!

    num_vertices++;
    if (num_vertices > vertices.rows())
      vertices.conservativeResize(vertices.rows()+10000, Eigen::NoChange);

    vertices.row(num_vertices-1)  = midpt;
    edge2vertex[EdgeKey1(i0, i1)] = num_vertices-1;

    return num_vertices-1;
  }
  ;

  // maps an edge to the sample vertex generated on it
  MyMap  edge2vertex;
};


IGL_INLINE void igl::copyleft::marching_cubes_root_finding(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const double R,
  const double cutoff,
  const Eigen::MatrixXd &P,
  const Eigen::MatrixXd &points,
  const unsigned x_res,
  const unsigned y_res,
  const unsigned z_res,
  const igl::copyleft::LocalImplicitFunction LOCAL_IMPLICIT_FUNCTION_TYPE,
  Eigen::MatrixXd &vertices,
  Eigen::MatrixXi &faces)
{
  MarchingCubesRF mc(V, F, R, cutoff, P,
                                       points,
                                       x_res,
                                       y_res,
                                       z_res,
                                       LOCAL_IMPLICIT_FUNCTION_TYPE,
                                       vertices,
                                       faces);
}
