#include "derive.hpp"

#include <iostream>

#include "access.hpp"
#include "loop.hpp"
#include "Omega_h_map.hpp"
#include "project.hpp"
#include "space.hpp"

namespace Omega_h {

template <Int dim>
INLINE Matrix<dim, dim> get_simplex_jacobian(Few<Vector<dim>, dim + 1> evv2x) {
  Matrix<dim, dim> dx_dxi;
  for (Int i = 0; i < dim; ++i) {
    for (Int j = 0; j < dim; ++j) {
      dx_dxi[j][i] = evv2x[i + 1][j] - evv2x[0][j];
    }
  }
  return dx_dxi;
}

template <Int dim>
static Reals derive_element_gradients_dim(Mesh* mesh, Reals vert_values) {
  auto coords = mesh->coords();
  auto ev2v = mesh->ask_elem_verts();
  auto out = Write<Real>(mesh->nelems() * dim);
  auto f = LAMBDA(LO e) {
    auto evv2v = gather_verts<dim + 1>(ev2v, e);
    auto evv2u = gather_scalars<dim + 1>(vert_values, evv2v);
    Vector<dim> du_dxi;
    for (Int i = 0; i < dim; ++i) du_dxi[i] = evv2u[i + 1] - evv2u[0];
    auto evv2x = gather_vectors<dim + 1, dim>(coords, evv2v);
    Matrix<dim, dim> dx_dxi = get_simplex_jacobian<dim>(evv2x);
    auto dxi_dx = invert(dx_dxi);
    auto du_dx = dxi_dx * du_dxi;
    set_vector(out, e, du_dx);
  };
  parallel_for(mesh->nelems(), f);
  return out;
}

template <Int dim>
static Reals derive_element_hessians_dim(Mesh* mesh, Reals vert_gradients) {
  auto coords = mesh->coords();
  auto ev2v = mesh->ask_elem_verts();
  auto out = Write<Real>(mesh->nelems() * symm_dofs(dim));
  auto f = LAMBDA(LO e) {
    auto evv2v = gather_verts<dim + 1>(ev2v, e);
    auto evv2u = gather_vectors<dim + 1, dim>(vert_gradients, evv2v);
    Matrix<dim, dim> du_dxi;
    for (Int i = 0; i < dim; ++i) {
      for (Int j = 0; j < dim; ++j) {
        du_dxi[i][j] = evv2u[j + 1][i] - evv2u[0][i];
      }
    }
    auto evv2x = gather_vectors<dim + 1, dim>(coords, evv2v);
    auto dx_dxi = get_simplex_jacobian<dim>(evv2x);
    auto dxi_dx = invert(dx_dxi);
    auto du_dx = dxi_dx * du_dxi;
    set_symm(out, e, du_dx);
  };
  parallel_for(mesh->nelems(), f);
  return out;
}

Reals derive_element_gradients(Mesh* mesh, Reals vert_values) {
  CHECK(vert_values.size() == mesh->nverts());
  if (mesh->dim() == 3)
    return derive_element_gradients_dim<3>(mesh, vert_values);
  if (mesh->dim() == 2)
    return derive_element_gradients_dim<2>(mesh, vert_values);
  NORETURN(Reals());
}

Reals derive_element_hessians(Mesh* mesh, Reals vert_gradients) {
  CHECK(vert_gradients.size() == mesh->nverts() * mesh->dim());
  if (mesh->dim() == 3)
    return derive_element_hessians_dim<3>(mesh, vert_gradients);
  if (mesh->dim() == 2)
    return derive_element_hessians_dim<2>(mesh, vert_gradients);
  NORETURN(Reals());
}

Reals recover_gradients(Mesh* mesh, Reals vert_values) {
  auto e_grad = derive_element_gradients(mesh, vert_values);
  return project_by_fit(mesh, e_grad);
}

Reals recover_hessians_from_gradients(Mesh* mesh, Reals vert_gradients) {
  auto e_hess = derive_element_hessians(mesh, vert_gradients);
  return project_by_fit(mesh, e_hess);
}

Reals recover_hessians(Mesh* mesh, Reals vert_values) {
  auto v_grad = recover_gradients(mesh, vert_values);
  return recover_hessians_from_gradients(mesh, v_grad);
}
}  // namespace Omega_h
