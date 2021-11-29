/// \file half_edge_mesh_coord.hpp
/// template classes for coordinate computations 
///   (midpoints, edge lengths, angles, etc.)

/*
  Copyright (C) 2021 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <limits>
#include <cmath>

// *** DEBUG ***
// Uncomment if you want to include print statements 
//   for debugging in this file.
#include<iostream>


#ifndef _HALF_EDGE_MESH_COORD_
#define _HALF_EDGE_MESH_COORD_

namespace HMESH {

  /// Compute the midpoint of coord0[] and coord1[] in coord2[].
  template <const int DIM, typename CTYPE>
  void compute_midpoint
  (const CTYPE coord0[DIM], const CTYPE coord1[DIM],
   CTYPE coord2[DIM])
  {
    for (int d = 0; d < DIM; d++)
      { coord2[d] = (coord0[d] + coord1[d])/2.0; };
  }


  /// Compute the squared distance between two coordinates.
  template <const int DIM, typename CTYPE>
  double compute_squared_distance
  (const CTYPE coord0[DIM], const CTYPE coord1[DIM])
  {
    double sum = 0.0;
    for (int d = 0; d < DIM; d++) {
      const double diff = coord0[d] - coord1[d];
      sum += (diff*diff);
    }

    return sum;
  }


  /// Compute magnitude squared.
  template <const int DIM, typename CTYPE>
  double compute_magnitude_squared(const CTYPE coord[DIM])
  {
    double magnitude_squared = 0.0;
    for (int d = 0; d < DIM; d++) {
      CTYPE c = coord[d];
      magnitude_squared += (c*c);
    }

    return magnitude_squared;
  }


  /// Set all coordinates to c.
  template <const int DIM, typename CTYPEA, typename CTYPEB>
  void set_all_coord(const CTYPEA c, CTYPEB coord[DIM])
  {
    for (int d = 0; d < DIM; d++)
      { coord[d] = CTYPEB(c); }
  }


  /// Copy coord0 to coord1.
  template <const int DIM, typename CTYPE0, typename CTYPE1>
  void copy_coord(const CTYPE0 coord0[DIM], CTYPE1 coord1[DIM])
  {
    for (int d = 0; d < DIM; d++)
      { coord1[d] = CTYPE1(coord0[d]); }
  }
  
  /// Divide coord by a (non-zero) scalar.
  template <const int DIM, typename CTYPEA, typename CTYPEB>
  void divide_by_scalar(const CTYPEA c, CTYPEB coord[DIM])
  {
    for (int d = 0; d < DIM; d++)
      { coord[d] = coord[d]/c; }
  }


  /// Subtract coord1 from coord0[].
  /// @param[out] coord2[] = coord0[] - coord1[].
  template <const int DIM, typename CTYPE0, 
            typename CTYPE1, typename CTYPE2>
  void subtract_coord(const CTYPE0 coord0[DIM], CTYPE1 coord1[DIM],
                      CTYPE2 coord2[DIM])
  {
    for (int d = 0; d < DIM; d++)
      { coord2[d] = coord0[d] - coord1[d]; }
  }


  /// Normalize vector.
  /// - Also returns the magnitude of the original vector.
  /// - If all coordinates are 0, the magnitude will be
  ///   set to 0, and the vector to (1,0,0,...).
  /// @tparam DIM Dimension. Must be at least 1.
  template <const int DIM, typename CTYPE>
  void normalize_vector(CTYPE coord[DIM], double & magnitude)
  {
    const double magnitude_squared = 
      compute_magnitude_squared<DIM>(coord);

    magnitude = std::sqrt(magnitude_squared);

    if (std::abs(magnitude) == 0.0) {
      magnitude = 0.0;
      set_all_coord<DIM>(0, coord);
      coord[0] = 1.0;
      return;
    }

    double max_abs_coord = std::abs(coord[0]);
    int dmax = 0;
    
    // Extra processing to avoid overflow if vector has small magnitude.

    // Get the coordinate with maximum absolute value.
    for (int d = 1; d < DIM; d++) {
      double abs_c = std::abs(coord[d]);
      if (abs_c > max_abs_coord) {
        max_abs_coord = abs_c;
        dmax = d;
      }
    }

    if (max_abs_coord == 0) {
      magnitude = 0.0;
      set_all_coord<DIM>(0, coord);
      coord[0] = 1.0;
      return;
    }

    // Divide by max_abs_coord.
    for (int d = 0; d < DIM; d++) {
      if (d == dmax) {
        // Ensure that coord[d] is 1 or -1.
        if (coord[d] < 0) 
          { coord[d] = -1; }
        else 
          { coord[d] = 1; }
      }
      else {
        // Note: abs(coord[d]) <= max_abs_coord.
        coord[d] = coord[d]/max_abs_coord;
      }
    }

    // Since coord[dmax] is 1 or -1, magnitudeB >= 1.
    double magnitudeB = compute_magnitude_squared<DIM>(coord);
    magnitudeB = std::sqrt(magnitudeB);

    // Divide by magnitudeB.
    divide_by_scalar<DIM>(magnitudeB, coord);
  }


  /// Normalize vector.
  /// - Also returns the magnitude of the original vector.
  /// - Version that returns normalized vector in coord1[].
  template <const int DIM, typename CTYPE0, typename CTYPE1>
  void normalize_vector
  (const CTYPE0 coord0[DIM], CTYPE1 coord1[DIM], double & magnitude)
  {
    copy_coord<DIM>(coord0, coord1);
    normalize_vector<DIM>(coord1, magnitude);
  }


  /// Return the inner product (dot product) of two vectors.
  template <const int DIM, typename CTYPE>
  double compute_inner_product
  (const CTYPE coord0[DIM], const CTYPE coord1[DIM])
  {
    double product = 0.0;
    for (int d = 0; d < DIM; d++) 
      { product += (coord0[d] * coord1[d]); }

    return product;
  }

  /// Compute cosine of angle between two vectors.
  /// - If either vector is 0 (or very, very near 0),
  ///   sets flag_zero to true and returns 0.
  template <const int DIM, typename CTYPE>
  double compute_cos_angle
  (const CTYPE coord0[DIM], const CTYPE coord1[DIM], 
   bool & flag_zero)
  {
    double magnitude0, magnitude1;
    CTYPE temp_coord0[DIM];
    CTYPE temp_coord1[DIM];

    // Initialize.
    flag_zero = false;
    
    normalize_vector<DIM>(coord0, temp_coord0, magnitude0);
    normalize_vector<DIM>(coord1, temp_coord1, magnitude1);

    if ((magnitude0 == 0.0) || (magnitude1 == 0.0)) {
      flag_zero = true;
      return 0;
    }

    double cos_angle =
      compute_inner_product<DIM>(temp_coord0, temp_coord1);

    // Clamp to [-1,1] to handle numerical errors.
    if (cos_angle < -1) { cos_angle = -1; }
    if (cos_angle > 1) { cos_angle = 1; }

    return cos_angle;
  }


  /// Compute cos of triangle angle at coord1[] 
  ///   in triangle (coord0[], coord1[], coord2[]).
  /// - If coord0[] == coord1[] (or is very, very close,) or
  ///   coord2[] == coord1[] (or is very, very close,)
  ///   sets flag_zero to true and returns 0.
  template <const int DIM, typename CTYPE>
  double compute_cos_triangle_angle
  (const CTYPE coord0[DIM], const CTYPE coord1[DIM], const CTYPE coord2[],
   bool & flag_zero) 
  {
    CTYPE vectA[DIM];
    CTYPE vectB[DIM];

    subtract_coord<DIM>(coord0, coord1, vectA);
    subtract_coord<DIM>(coord2, coord1, vectB);

    return compute_cos_angle<DIM>(vectA, vectB, flag_zero);
  }
  
}

#endif
