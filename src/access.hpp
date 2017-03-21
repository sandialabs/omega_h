#ifndef ACCESS_HPP
#define ACCESS_HPP

#include "Omega_h_math.hpp"
#include "algebra.hpp"
#include "internal.hpp"

namespace Omega_h {

template <Int nhhl>
DEVICE Few<LO, nhhl> gather_down(LOs const& hl2l, Int h) {
  Few<LO, nhhl> hhl2l;
  for (Int i = 0; i < nhhl; ++i) {
    auto hl = h * nhhl + i;
    hhl2l[i] = hl2l[hl];
  }
  return hhl2l;
}

template <Int neev>
DEVICE Few<LO, neev> gather_verts(LOs const& ev2v, Int e) {
  return gather_down<neev>(ev2v, e);
}

template <Int neev, typename T>
DEVICE Few<T, neev> gather_scalars(Read<T> const& a, Few<LO, neev> v) {
  Few<T, neev> x;
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
  set_vector(a, i, matrix2vector(m));
}

template <Int dim>
DEVICE Matrix<dim, dim> get_matrix(Reals const& a, Int i) {
  return vector2symm(get_vector<matrix_dofs(dim)>(a, i));
}

Reals vectors_2d_to_3d(Reals vecs2);
Reals vectors_3d_to_2d(Reals vecs2);

Reals average_field(Mesh* mesh, Int dim, LOs a2e, Int ncomps, Reals v2x);

}  // end namespace Omega_h

#endif
