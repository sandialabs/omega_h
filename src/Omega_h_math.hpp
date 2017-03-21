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

template <Int n>
OMEGA_H_DEVICE void set_vector(Write<Real> const& a, Int i, Vector<n> v) {
  for (Int j = 0; j < n; ++j) a[i * n + j] = v[j];
}

template <Int n, class Arr>
OMEGA_H_DEVICE Vector<n> get_vector(Arr const& a, Int i) {
  Vector<n> v;
  for (Int j = 0; j < n; ++j) v[j] = a[i * n + j];
  return v;
}

template <Int n>
OMEGA_H_DEVICE void set_symm(Write<Real> const& a, Int i, Matrix<n, n> symm) {
  set_vector(a, i, symm2vector(symm));
}

template <Int n, typename Arr>
OMEGA_H_DEVICE Matrix<n, n> get_symm(Arr const& a, Int i) {
  return vector2symm(get_vector<symm_dofs(n)>(a, i));
}

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
