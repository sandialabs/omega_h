#include "Omega_h_recover.hpp"

/* This code projects an element-wise
 * field onto the vertices.
 * All interior vertices fit a linear polynomial
 * based on surrounding element values
 * and evaluate that polynomial to obtain their
 * values.
 * Boundary nodes evaluate the polynomial of a
 * nearby interior node. This involves storing
 * the polynomial coefficients on the interior nodes
 * and copying them over to the boundary in a diffusive
 * fashion.
 *
 * This approach is inspired by the following paper:
 * Boussetta, Ramzy, Thierry Coupez, and Lionel Fourment.
 * "Adaptive remeshing based on a posteriori error estimation
 *  for forging simulation."
 * Computer methods in applied mechanics and engineering
 * 195.48 (2006): 6626-6645.
 */

#include "Omega_h_array_ops.hpp"
#include "Omega_h_fit.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

template <Int dim>
Reals get_interior_coeffs_dim(Mesh* mesh, Reals e_data, Int ncomps) {
  auto v2e = mesh->ask_up(VERT, dim);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ev2v = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto owned = mesh->owned(VERT);
  auto class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto out = Write<Real>(mesh->nverts() * ncomps * (dim + 1));
  auto f = OMEGA_H_LAMBDA(LO v) {
    if (!owned[v] || (class_dim[v] != dim)) return;
    auto qr = get_cavity_qr_factorization<dim>(v, v2ve, ve2e, ev2v, coords);
    for (Int comp = 0; comp < ncomps; ++comp) {
      auto coeffs =
          fit_cavity_polynomial<dim>(qr, v, v2ve, ve2e, e_data, comp, ncomps);
      set_vector(out, v * ncomps + comp, coeffs);
    }
  };
  parallel_for(mesh->nverts(), f, "get_interior_coeffs");
  return mesh->sync_array(VERT, Reals(out), ncomps * (dim + 1));
}

static Reals get_interior_coeffs(Mesh* mesh, Reals e_data, Int ncomps) {
  if (mesh->dim() == 3) {
    return get_interior_coeffs_dim<3>(mesh, e_data, ncomps);
  } else if (mesh->dim() == 2) {
    return get_interior_coeffs_dim<2>(mesh, e_data, ncomps);
  } else if (mesh->dim() == 1) {
    return get_interior_coeffs_dim<1>(mesh, e_data, ncomps);
  }
  OMEGA_H_NORETURN(Reals());
}

void diffuse_to_exterior(
    Mesh* mesh, Reals* p_v_data, Int ncomps, Read<I8>* p_visited) {
  auto v_data = *p_v_data;
  auto visited = *p_visited;
  auto new_data = deep_copy(v_data);
  auto new_visited = deep_copy(visited);
  auto v2v = mesh->ask_star(VERT);
  auto v2vv = v2v.a2ab;
  auto vv2v = v2v.ab2b;
  auto f = OMEGA_H_LAMBDA(LO v) {
    if (visited[v]) return;
    Int nadj = 0;
    for (auto vv = v2vv[v]; vv < v2vv[v + 1]; ++vv) {
      auto ov = vv2v[vv];
      if (visited[ov]) ++nadj;
    }
    if (!nadj) return;
    for (Int comp = 0; comp < ncomps; ++comp) {
      Real sum = 0.0;
      for (auto vv = v2vv[v]; vv < v2vv[v + 1]; ++vv) {
        auto ov = vv2v[vv];
        if (visited[ov]) sum += v_data[ov * ncomps + comp];
      }
      new_data[v * ncomps + comp] = sum / nadj;
    }
    new_visited[v] = 1;
  };
  parallel_for(mesh->nverts(), f, "diffuse_to_exterior");
  v_data = new_data;
  visited = new_visited;
  v_data = mesh->sync_array(VERT, v_data, ncomps);
  visited = mesh->sync_array(VERT, visited, 1);
  *p_v_data = v_data;
  *p_visited = visited;
}

template <Int dim>
Reals evaluate_coeffs_dim(Mesh* mesh, Reals v_coeffs, Int ncomps) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * ncomps);
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto x = get_vector<dim>(coords, v);
    for (Int comp = 0; comp < ncomps; ++comp) {
      auto coeffs = get_vector<dim + 1>(v_coeffs, v * ncomps + comp);
      auto val = eval_polynomial(coeffs, x);
      out[v * ncomps + comp] = val;
    }
  };
  parallel_for(mesh->nverts(), f, "evaluate_coeffs");
  return out;
}

static Reals evaluate_coeffs(Mesh* mesh, Reals v_coeffs, Int ncomps) {
  if (mesh->dim() == 3) {
    return evaluate_coeffs_dim<3>(mesh, v_coeffs, ncomps);
  } else if (mesh->dim() == 2) {
    return evaluate_coeffs_dim<2>(mesh, v_coeffs, ncomps);
  } else if (mesh->dim() == 1) {
    return evaluate_coeffs_dim<1>(mesh, v_coeffs, ncomps);
  }
  OMEGA_H_NORETURN(Reals());
}

bool has_interior_verts(Mesh* mesh) {
  auto dim = mesh->dim();
  auto class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto interior = each_eq_to(class_dim, I8(dim));
  auto have_local_interior = (get_max(interior) == 1);
  auto comm = mesh->comm();
  return comm->reduce_or(have_local_interior);
}

Reals project_by_fit(Mesh* mesh, Reals e_data) {
  OMEGA_H_CHECK(mesh->owners_have_all_upward(VERT));
  OMEGA_H_CHECK(e_data.size() % mesh->nelems() == 0);
  OMEGA_H_CHECK(has_interior_verts(mesh));
  auto ncomps = e_data.size() / mesh->nelems();
  auto dim = mesh->dim();
  auto v_coeffs = get_interior_coeffs(mesh, e_data, ncomps);
  auto class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto visited = each_eq_to(class_dim, I8(dim));
  while (mesh->comm()->reduce_or(get_min(visited) == 0)) {
    diffuse_to_exterior(mesh, &v_coeffs, ncomps * (dim + 1), &visited);
  }
  return evaluate_coeffs(mesh, v_coeffs, ncomps);
}

Reals project_by_average(Mesh* mesh, Reals e_data) {
  OMEGA_H_CHECK(mesh->owners_have_all_upward(VERT));
  OMEGA_H_CHECK(e_data.size() % mesh->nelems() == 0);
  auto ncomps = e_data.size() / mesh->nelems();
  auto verts2elems = mesh->ask_up(VERT, mesh->dim());
  auto weights = Reals(verts2elems.ab2b.size(), 1.0);
  auto avgs = graph_weighted_average(verts2elems, weights, e_data, ncomps);
  return mesh->sync_array(VERT, avgs, ncomps);
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> get_simplex_jacobian(
    Few<Vector<dim>, dim + 1> evv2x) {
  Matrix<dim, dim> dx_dxi;
  for (Int i = 0; i < dim; ++i) {
    for (Int j = 0; j < dim; ++j) {
      dx_dxi[j][i] = evv2x[i + 1][j] - evv2x[0][j];
    }
  }
  return dx_dxi;
}

template <Int dim>
Reals derive_element_gradients_dim(Mesh* mesh, Reals vert_values) {
  auto coords = mesh->coords();
  auto ev2v = mesh->ask_elem_verts();
  auto out = Write<Real>(mesh->nelems() * dim);
  auto f = OMEGA_H_LAMBDA(LO e) {
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
  parallel_for(mesh->nelems(), f, "derive_element_gradients");
  return out;
}

template <Int dim>
Reals derive_element_hessians_dim(Mesh* mesh, Reals vert_gradients) {
  auto coords = mesh->coords();
  auto ev2v = mesh->ask_elem_verts();
  auto out = Write<Real>(mesh->nelems() * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO e) {
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
  parallel_for(mesh->nelems(), f, "derive_element_hessians");
  return out;
}

Reals derive_element_gradients(Mesh* mesh, Reals vert_values) {
  OMEGA_H_CHECK(vert_values.size() == mesh->nverts());
  if (mesh->dim() == 3) {
    return derive_element_gradients_dim<3>(mesh, vert_values);
  } else if (mesh->dim() == 2) {
    return derive_element_gradients_dim<2>(mesh, vert_values);
  } else if (mesh->dim() == 1) {
    return derive_element_gradients_dim<1>(mesh, vert_values);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals derive_element_hessians(Mesh* mesh, Reals vert_gradients) {
  OMEGA_H_CHECK(vert_gradients.size() == mesh->nverts() * mesh->dim());
  if (mesh->dim() == 3) {
    return derive_element_hessians_dim<3>(mesh, vert_gradients);
  } else if (mesh->dim() == 2) {
    return derive_element_hessians_dim<2>(mesh, vert_gradients);
  } else if (mesh->dim() == 1) {
    return derive_element_hessians_dim<1>(mesh, vert_gradients);
  }
  OMEGA_H_NORETURN(Reals());
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

}  // end namespace Omega_h
