#include <Omega_h_assemble.hpp>

#include <Omega_h_align.hpp>
#include <Omega_h_element.hpp>
#include <Omega_h_loop.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>

namespace Omega_h {

template <Int dim>
Reals get_edge_grad_grad_dim(Mesh* mesh, Reals elem_jac_invs,
    Reals elem_material_matrices, Reals elem_sizes) {
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
      auto elem_vert0 =
          element_down_template(family, dim, EDGE, elem_edge, rot);
      auto elem_vert1 =
          element_down_template(family, dim, EDGE, elem_edge, 1 - rot);
      auto jac_inv = get_matrix<dim>(elem_jac_invs, elem);
      auto grad0 =
          (elem_vert0 == 0) ? (-sum(jac_inv)) : jac_inv[elem_vert0 - 1];
      auto grad1 =
          (elem_vert1 == 0) ? (-sum(jac_inv)) : jac_inv[elem_vert1 - 1];
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
Reals get_vert_grad_grad_dim(Mesh* mesh, Reals elem_jac_invs,
    Reals elem_material_matrices, Reals elem_sizes) {
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

Reals get_edge_grad_grad(Mesh* mesh, Reals elem_jac_invs,
    Reals elem_material_matrices, Reals elem_sizes) {
  if (mesh->dim() == 3)
    return get_edge_grad_grad_dim<3>(
        mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  if (mesh->dim() == 2)
    return get_edge_grad_grad_dim<2>(
        mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  if (mesh->dim() == 1)
    return get_edge_grad_grad_dim<1>(
        mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  OMEGA_H_NORETURN(Reals());
}

Reals get_vert_grad_grad(Mesh* mesh, Reals elem_jac_invs,
    Reals elem_material_matrices, Reals elem_sizes) {
  if (mesh->dim() == 3)
    return get_vert_grad_grad_dim<3>(
        mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  if (mesh->dim() == 2)
    return get_vert_grad_grad_dim<2>(
        mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  if (mesh->dim() == 1)
    return get_vert_grad_grad_dim<1>(
        mesh, elem_jac_invs, elem_material_matrices, elem_sizes);
  OMEGA_H_NORETURN(Reals());
}

Reals get_elem_jac_invs(Mesh* mesh) {
  if (mesh->dim() == 3) return get_elem_jac_invs_dim<3>(mesh);
  if (mesh->dim() == 2) return get_elem_jac_invs_dim<2>(mesh);
  if (mesh->dim() == 1) return get_elem_jac_invs_dim<1>(mesh);
  OMEGA_H_NORETURN(Reals());
}

/*
   A*x = b
   B = A*e_i
   (A - A*e_i)*x = A*x - A*e_i*x
   (A - A*e_i)*x = b - A*e_i*x
 */

void apply_dirichlet(Mesh* mesh, Reals* p_a_edge, Reals a_vert, Reals* p_b,
    Bytes are_dbc, Reals dbc_values) {
  auto a_edge = *p_a_edge;
  auto b = *p_b;
  OMEGA_H_CHECK(a_edge.size() == mesh->nedges());
  OMEGA_H_CHECK(a_vert.size() == mesh->nverts());
  auto a_edge_new = deep_copy(a_edge);
  auto b_new = deep_copy(b);
  auto verts2edges = mesh->ask_up(VERT, EDGE);
  auto edges2verts = mesh->ask_verts_of(EDGE);
  auto v_f = OMEGA_H_LAMBDA(LO v) {
    Real b_v_new = b_new[v];
    if (are_dbc[v]) {
      b_v_new = a_vert[v] * dbc_values[v];
    } else {
      auto begin = verts2edges.a2ab[v];
      auto end = verts2edges.a2ab[v + 1];
      for (auto vert_edge = begin; vert_edge < end; ++vert_edge) {
        auto edge = verts2edges.ab2b[vert_edge];
        auto code = verts2edges.codes[vert_edge];
        auto edge_vert = code_which_down(code);
        auto other_vert = edges2verts[edge * 2 + (1 - edge_vert)];
        if (are_dbc[other_vert]) {
          b_v_new -= a_edge[edge] * dbc_values[other_vert];
        }
      }
    }
    b_new[v] = b_v_new;
  };
  parallel_for(mesh->nverts(), v_f);
  *p_b = mesh->sync_array(VERT, Reals(b_new), 1);
  auto e_f = OMEGA_H_LAMBDA(LO e) {
    auto a_e_new = a_edge_new[e];
    for (Int edge_vert = 0; edge_vert < 2; ++edge_vert) {
      if (are_dbc[edges2verts[e * 2 + edge_vert]]) {
        a_e_new = 0.;
      }
    }
    a_edge_new[e] = a_e_new;
  };
  parallel_for(mesh->nedges(), e_f);
  *p_a_edge = mesh->sync_array(EDGE, Reals(a_edge_new), 1);
}

}  // namespace Omega_h
