// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_MARCHINGCUBESROOTFINDING_H
#define IGL_COPYLEFT_MARCHINGCUBESROOTFINDING_H
#include "../igl_inline.h"

#include <Eigen/Core>
namespace igl
{
  namespace copyleft
  {
    // marching_cubes_root_finding(R, P, points, x_res, y_res, z_res, vertices, faces )
    //
    // performs marching cubes reconstruction by root finding on a implicit function parametrized 
    // by the point cloud
    //
    // Note the current tol for root finding is fixed at 1e-4 for the sqaured norm distance between two subsequent approximations of the root.
    //
    //
    // Citation: Wyvill et al., "Data structure for soft objects", The Visual Computer, August 1986, Volume 2, Issue 4, pp 227–234.
    //           http://thesesups.ups-tlse.fr/3962/1/2018TOU30055.pdf
    //
    // Input:
    //  R Radius of influence
    //  cutoff Value of the field for points on a surface. If < 0 then function value at each of the input point cloud P is calculated 
    //         the radius of convergence and averaged to get the cutoff (similar to Mesh Reconstruction HW).
    //  P  #P by 3 list of input points of the point cloud
    //  points  #number_of_grid_points x 3 array -- 3-D positions of the grid
    //    points, ordered in x,y,z order:
    //      points[index] = the point at (x,y,z) where :
    //      x = (index % (xres -1),
    //      y = (index / (xres-1)) %(yres-1),
    //      z = index / (xres -1) / (yres -1) ).
    //      where x,y,z index x, y, z dimensions
    //      i.e. index = x + y*xres + z*xres*yres
    //  xres  resolutions of the grid in x dimension
    //  yres  resolutions of the grid in y dimension
    //  zres  resolutions of the grid in z dimension
    // Output:
    //   vertices  #V by 3 list of mesh vertex positions
    //   faces  #F by 3 list of mesh triangle indices
    //
    enum LocalImplicitFunction
    {
      LOCAL_IMPLICIT_FUNCTION_DEFAULT = 0, // same as WYVILL_1986
      LOCAL_IMPLICIT_FUNCTION_WYVILL_1986 = 1
    };
    
    template <
      typename Derivedvalues, 
      typename Derivedpoints, 
      typename Derivedvertices, 
      typename DerivedF>
      IGL_INLINE void marching_cubes_root_finding(
        const double R,
        const double cutoff,
        const Eigen::PlainObjectBase<Derivedvalues> &P,
        const Eigen::PlainObjectBase<Derivedpoints> &points,
        const unsigned x_res,
        const unsigned y_res,
        const unsigned z_res,
        const LocalImplicitFunction LOCAL_IMPLICIT_FUNCTION_TYPE,
        Eigen::PlainObjectBase<Derivedvertices> &vertices,
        Eigen::PlainObjectBase<DerivedF> &faces);
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "marching_cubes_root_finding.cpp"
#endif

#endif
