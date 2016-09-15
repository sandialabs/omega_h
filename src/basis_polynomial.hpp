#ifndef BASIS_POLYNOMIAL_HPP
#define BASIS_POLYNOMIAL_HPP

#include "Omega_h_r3d.hpp"
#include "size.hpp"

namespace Omega_h {

/* Convert one of the linear simplex basis
 * functions (vertex basis functions) into
 * a linear polynomial in global (x,y,z) coordinates.
 * Here we distinguish between (xi) which is the
 * vector of 3 independent barycentric coordinates,
 * and (eta), which represents one of the four basis functions.
 * as the (dxi_deta) logic suggests, (eta) may be equal
 * to one of the three (xi) components or it is the
 * leftover barycentric coordinate defined by
 *   eta = 1 - xi[0] - xi[1] - xi[2];
 */

template <Int dim>
INLINE r3d::Polynomial<dim, 1> get_basis_polynomial(
    Few<Vector<dim>, dim + 1> elem_pts,
    Int elem_vert) {
  auto dx_dxi = simplex_basis<dim, dim>(elem_pts);
  auto dxi_dx = invert(J);
  Vector<3> dxi_deta;
  if (elem_vert) {
    dxi_deta = zero_vector<dim>();
    dxi_deta[elem_vert - 1] = 1;
  } else {
    for (Int i = 0; i < dim; ++i)
      dxi_deta[i] = -1;
  }
  auto dx_deta = dx_dxi * dxi_deta;
  auto deta_dx = pseudo_invert(dx_deta);
  auto other_vert = (elem_vert + 1) % (dim + 1);
  auto origin = elem_pts[other_vert];
  auto origin_val = deta_dx * (-origin);
  r3d::Polynomial<dim, 1> poly;
  poly.coeffs[0] = origin_val;
  for (Int i = 0; i < dim; ++i)
    poly.coeffs[i + 1] = deta_dx[i];
  return poly;
}

}

#endif
