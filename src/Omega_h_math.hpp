#ifndef OMEGA_H_MATH_HPP
#define OMEGA_H_MATH_HPP

#include <cfloat>
#include <climits>
#include <cmath>

#include <Omega_h_array.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_few.hpp>
#include <Omega_h_scalar.hpp>
#include <Omega_h_vector.hpp>
#include <Omega_h_matrix.hpp>

namespace Omega_h {

template <Int dim>
OMEGA_H_INLINE Vector<dim> metric_eigenvalues(Vector<dim> h) {
  Vector<dim> l;
  for (Int i = 0; i < dim; ++i) l[i] = 1.0 / square(h[i]);
  return l;
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> compose_metric(
    Matrix<dim, dim> r, Vector<dim> h) {
  auto l = metric_eigenvalues(h);
  return r * diagonal(l) * transpose(r);
}

OMEGA_H_INLINE Real metric_eigenvalue_from_length(Real h) {
  return 1.0 / square(h);
}

}  // end namespace Omega_h

#endif
