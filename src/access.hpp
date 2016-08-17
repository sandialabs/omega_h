#ifndef ACCESS_HPP
#define ACCESS_HPP

#include "algebra.hpp"
#include "internal.hpp"

namespace Omega_h {

template <Int n, class Arr>
DEVICE Vector<n> get_vector(Arr const& a, Int i) {
  Vector<n> v;
  for (Int j = 0; j < n; ++j) v[j] = a[i * n + j];
  return v;
}

INLINE Matrix<2, 2> vector2symm(Vector<3> v) {
  Matrix<2, 2> symm;
  symm[0][0] = v[0];
  symm[1][1] = v[1];
  symm[1][0] = v[2];
  symm[0][1] = symm[1][0];
  return symm;
}

INLINE Matrix<3, 3> vector2symm(Vector<6> v) {
  Matrix<3, 3> symm;
  symm[0][0] = v[0];
  symm[1][1] = v[1];
  symm[2][2] = v[2];
  symm[1][0] = v[3];
  symm[2][1] = v[4];
  symm[2][0] = v[5];
  symm[0][1] = symm[1][0];
  symm[1][2] = symm[2][1];
  symm[0][2] = symm[2][0];
  return symm;
}

template <Int n, typename Arr>
DEVICE Matrix<n, n> get_symm(Arr const& a, Int i) {
  return vector2symm(get_vector<symm_dofs(n)>(a, i));
}

template <Int neev>
DEVICE Few<LO, neev> gather_verts(LOs const& ev2v, Int e) {
  Few<LO, neev> v;
  for (Int i = 0; i < neev; ++i) v[i] = ev2v[e * neev + i];
  return v;
}

template <Int neev>
DEVICE Few<Real, neev> gather_scalars(Reals const& a, Few<LO, neev> v) {
  Few<Real, neev> x;
  for (Int i = 0; i < neev; ++i) x[i] = a[v[i]];
  return x;
}

template <Int neev, Int dim>
DEVICE Few<Vector<dim>, neev> gather_vectors(Reals const& a, Few<LO, neev> v) {
  Few<Vector<dim>, neev> x;
  for (Int i = 0; i < neev; ++i) x[i] = get_vector<dim>(a, v[i]);
  return x;
}

template <Int neev, Int dim>
DEVICE Few<Matrix<dim, dim>, neev> gather_symms(
    Reals const& a, Few<LO, neev> v) {
  Few<Matrix<dim, dim>, neev> x;
  for (Int i = 0; i < neev; ++i) x[i] = get_symm<dim>(a, v[i]);
  return x;
}

template <Int dim>
DEVICE void set_matrix(Write<Real> const& a, Int i, Matrix<dim, dim> m) {
  for (Int j = 0; j < dim; ++j) {
    for (Int k = 0; k < dim; ++k) {
      a[(i * dim + j) * dim + k] = m[j][k];
    }
  }
}

template <Int dim>
DEVICE Matrix<dim, dim> get_matrix(Reals const& a, Int i) {
  Matrix<dim, dim> m;
  for (Int j = 0; j < dim; ++j) {
    for (Int k = 0; k < dim; ++k) {
      m[j][k] = a[(i * dim + j) * dim + k];
    }
  }
  return m;
}

Reals vectors_2d_to_3d(Reals vecs2);
Reals vectors_3d_to_2d(Reals vecs2);

Reals average_field(Mesh* mesh, Int dim, LOs a2e, Int ncomps, Reals v2x);

}  // end namespace Omega_h

#endif
