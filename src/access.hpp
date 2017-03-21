#ifndef ACCESS_HPP
#define ACCESS_HPP

#include "Omega_h_math.hpp"
#include "algebra.hpp"
#include "internal.hpp"

namespace Omega_h {

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
