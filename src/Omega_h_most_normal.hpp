#ifndef OMEGA_H_MOST_NORMAL_HPP
#define OMEGA_H_MOST_NORMAL_HPP

#include <Omega_h_vector.hpp>

namespace Omega_h {

/*
Aubry, Romain, and Rainald Lohner.
"On the 'most normal' normal"
Communications in Numerical Methods in Engineering
2008; 24:1641-1652.
DOI: 10.1002/cnm.1056
*/
template <Int nmax>
OMEGA_H_INLINE
Vector<3> get_most_normal_normal(Few<Vector<3>, nmax> N, Int n,
    Real inside_tol = OMEGA_H_EPSILON) {
  // first try pairs
  Real scalmin = -1.0;
  Vector<3> c;
  for (Int i = 0; i < n; ++i) {
    for (Int j = i + 1; j < n; ++j) {
      //Compute the bisector vector (center):
      auto N_b = N[i] + N[j];
      N_b = normalize(N_b);
      //Compute the radius (angle):
      auto scal = N_b * N[i];
      //Compare the current length against current
      //smallest length
      if (scal < scalmin) continue;
      //Are all the points inside this circle?
      Int l;
      for (l = 0; l < n; ++l) {
        auto scalt = N[l] * N_b;
        if (scalt + inside_tol < scal) break;
      }
      if (l < n) continue; // not all points inside
      //Store minimal radius and the center
      c = N_b;
      scalmin=scal;
    }
  }
  // if a circle with two points is found to enclose all the points...
  // then this is the smallest circumscribed circle
  if (scalmin != -1.0) return c;
  // no pair circle worked, try triplets
  for (Int i = 0; i < n; ++i) {
    for (Int j = i + 1; j < n; ++j) {
      for (Int k = j + 1; k < n; ++k) {
        //Compute the center of the circle
        auto r_ij = N[i] - N[j];
        auto rx_ij = r_ij[0];
        auto ry_ij = r_ij[1];
        auto rz_ij = r_ij[2];
        auto r_ik = N[i] - N[k];
        auto rx_ik = r_ik[0];
        auto ry_ik = r_ik[1];
        auto rz_ik = r_ik[2];
        auto denom = rx_ij * ry_ik - rx_ik * ry_ij;
        auto N_cx = (rz_ik * ry_ij - rz_ij * ry_ik) / denom;
        auto N_cy = (rx_ik * rz_ij - rx_ij * rz_ik) / denom;
        auto N_cz = 1.0 / std::sqrt(1.0 + N_cx * N_cx + N_cy * N_cy);
        N_cx = N_cx * N_cz;
        N_cy = N_cy * N_cz;
        auto N_c = vector_3(N_cx, N_cy, N_cz);
        //Check the orientation
        auto scal = N_c * N[i];
        if (scal < 0.0) N_c = - N_c;
        //Compute the radius
        scal = N_c * N[i];
        OMEGA_H_CHECK(are_close(scal, N_c * N[j]));
        OMEGA_H_CHECK(are_close(scal, N_c * N[k]));
        //Compare the current against smallest length
        if (scal < scalmin) continue;
        //Are all the points inside this circle?
        Int l;
        for (l = 0; l < n; ++l) {
          auto scalt = N[l] * N_c;
          if (scalt + inside_tol < scal) break;
        }
        if (l < n) continue; // not all points inside
        //Store the minimal radius and the center
        c = N_c;
        scalmin = scal;
      }
    }
  }
  OMEGA_H_CHECK(scalmin > -1.0);
  return c;
}

template <Int nmax>
OMEGA_H_INLINE
Vector<2> get_most_normal_normal(Few<Vector<2>, nmax> N, Int) {
  return normalize(N[0] + N[1]);
}

}  // end namespace Omega_h

#endif
