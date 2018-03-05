#include <Omega_h_assemble.hpp>

#include <Omega_h_mesh.hpp>
#include <Omega_h_align.hpp>
#include <Omega_h_element.hpp>
#include <Omega_h_shape.hpp>
#include <Omega_h_loop.hpp>

namespace Omega_h {

template <Int dim>
Reals get_edge_grad_grad_dim(Mesh* mesh, Reals elem_jac_invs, Reals elem_material_matrices, Reals elem_sizes) {
  OMEGA_H_CHECK(mesh->family() == OMEGA_H_SIMPLEX);
  OMEGA_H_CHECK(mesh->owners_have_all_upward(EDGE));
  auto edges2elems = mesh->ask_up(EDGE, dim);
  auto family = mesh->family();
  auto out_w = Write<Real>(mesh->nedges());
  auto f = OMEGA_H_LAMBDA(LO edge) {
    auto begin = edges2elems.a2ab[edge];
    auto end = edges2elems.a2ab[edge + 1];
    Real value = 0.;
    for (auto edge_elem = begin; edge_elem < end; ++edge_elem) {
      auto elem = edges2elems.ab2b[edge_elem];
      auto code = edges2elems.codes[edge_elem];
      auto elem_edge = code_which_down(code);
      auto rot = code_rotation(code);
      auto elem_vert0 = element_down_template(family, dim, EDGE, elem_edge, rot);
      auto elem_vert1 = element_down_template(family, dim, EDGE, elem_edge, 1 - rot);
      auto jac_inv = get_matrix<dim>(elem_jac_invs, elem);
      auto grad0 = (elem_vert0 == 0) ? (-sum(jac_inv)) : jac_inv[elem_vert0 - 1];
      auto grad1 = (elem_vert1 == 0) ? (-sum(jac_inv)) : jac_inv[elem_vert1 - 1];
      auto material_matrix = get_symm<dim>(elem_material_matrices, elem);
      auto elem_value = grad0 * (material_matrix * grad1);
      value += elem_value * elem_sizes[elem];
    }
    out_w[edge] = value;
  };
  parallel_for(mesh->nedges(), f);
  return mesh->sync_array(EDGE, Reals(out_w), 1);
}

template <Int dim>
Reals get_vert_grad_grad_dim(Mesh* mesh, Reals elem_jac_invs, Reals elem_material_matrices, Reals elem_sizes) {
  OMEGA_H_CHECK(mesh->family() == OMEGA_H_SIMPLEX);
  OMEGA_H_CHECK(mesh->owners_have_all_upward(VERT));
  auto verts2elems = mesh->ask_up(VERT, dim);
  auto out_w = Write<Real>(mesh->nverts());
  auto f = OMEGA_H_LAMBDA(LO vert) {
    auto begin = verts2elems.a2ab[vert];
    auto end = verts2elems.a2ab[vert + 1];
    Real value = 0.;
    for (auto vert_elem = begin; vert_elem < end; ++vert_elem) {
      auto elem = verts2elems.ab2b[vert_elem];
      auto code = verts2elems.codes[vert_elem];
      auto elem_vert = code_which_down(code);
      auto jac_inv = get_matrix<dim>(elem_jac_invs, elem);
      auto grad = (elem_vert == 0) ? (-sum(jac_inv)) : jac_inv[elem_vert - 1];
      auto material_matrix = get_symm<dim>(elem_material_matrices, elem);
      auto contrib = grad * (material_matrix * grad);
      value += contrib * elem_sizes[elem];
    }
    out_w[vert] = value;
  };
  parallel_for(mesh->nverts(), f);
  return mesh->sync_array(VERT, Reals(out_w), 1);
}

template <Int dim>
Reals get_elem_jac_invs_dim(Mesh* mesh) {
  OMEGA_H_CHECK(mesh->family() == OMEGA_H_SIMPLEX);
  constexpr auto nelem_verts = simplex_degree(dim, VERT);
  auto elems2verts = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto out_w = Write<Real>(mesh->nelems() * matrix_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO elem) {
    auto elem_verts = gather_verts<nelem_verts>(elems2verts, elem);
    auto elem_coords = gather_vectors<nelem_verts, dim>(coords, elem_verts);
    auto jac = simplex_basis<dim, dim>(elem_coords);
    auto jac_inv = transpose(invert(jac));
    set_matrix<dim>(out_w, elem, jac_inv);
  };
  parallel_for(mesh->nelems(), f);
  return out_w;
}

Reals get_edge_grad_grad(Mesh* mesh, Reals elem_jac_invs, Reals elem_material_matrices, Reals elem_sizes) {
  if (mesh->dim() == 3) return get_edge_grad_grad_dim<3>(mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  if (mesh->dim() == 2) return get_edge_grad_grad_dim<2>(mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  if (mesh->dim() == 1) return get_edge_grad_grad_dim<1>(mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  OMEGA_H_NORETURN(Reals());
}

Reals get_vert_grad_grad(Mesh* mesh, Reals elem_jac_invs, Reals elem_material_matrices, Reals elem_sizes) {
  if (mesh->dim() == 3) return get_vert_grad_grad_dim<3>(mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  if (mesh->dim() == 2) return get_vert_grad_grad_dim<2>(mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  if (mesh->dim() == 1) return get_vert_grad_grad_dim<1>(mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  OMEGA_H_NORETURN(Reals());
}

Reals get_elem_jac_invs(Mesh* mesh) {
  if (mesh->dim() == 3) return get_elem_jac_invs_dim<3>(mesh);
  if (mesh->dim() == 2) return get_elem_jac_invs_dim<2>(mesh);
  if (mesh->dim() == 1) return get_elem_jac_invs_dim<1>(mesh);
  OMEGA_H_NORETURN(Reals());
}

}
