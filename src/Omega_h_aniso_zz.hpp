#ifndef OMEGA_H_ANISO_ZZ_HPP
#define OMEGA_H_ANISO_ZZ_HPP

#include <Omega_h_math.hpp>
#include "algebra.hpp"

namespace Omega_h {

OMEGA_H_INLINE Matrix<2, 2> get_iso_tri_inv() {
  return invert(matrix_2x2(
      -sqrt(3.0) / 2.0, sqrt(3.0) / 2.0,
            -1.0 / 2.0,      -1.0 / 2.0));
}

OMEGA_H_INLINE Matrix<3, 3> get_iso_tet_inv() {
  return invert(matrix_3x3(
      -sqrt(2.0 / 3.0),  sqrt(2.0 / 3.0), 0.0,
      -sqrt(2.0) / 3.0, -sqrt(2.0) / 3.0, 2.0 * sqrt(2.0) / 3.0,
            -1.0 / 3.0,       -1.0 / 3.0, -1.0 / 3.0));
}

template <Int dim>
struct IsoSimplexInv;
template <>
struct IsoSimplexInv<2> {
  static OMEGA_H_INLINE Matrix<2, 2> get() { return get_iso_tri_inv(); }
};
template <>
struct IsoSimplexInv<3> {
  static OMEGA_H_INLINE Matrix<3, 3> get() { return get_iso_tet_inv(); }
};
template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> get_iso_simplex_inv() {
  return IsoSimplexInv<dim>::get();
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim>
get_centered_basis(Few<Vector<dim>, dim + 1> p) {
  auto center = average(p);
  Matrix<dim, dim> out;
  for (Int i = 0; i < dim; ++i) out[i] = p[i] - center;
  return out;
}

};

#endif
