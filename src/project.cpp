#include "project.hpp"

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

#include "array.hpp"
#include "fit.hpp"
#include "loop.hpp"

namespace osh {

template <Int dim>
static Reals get_interior_coeffs_dim(Mesh* mesh, Reals e_data, Int ncomps,
    Read<I8> interior) {
  auto v2e = mesh->ask_up(VERT, dim);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ev2v = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto owned = mesh->owned(VERT);
  auto out = Write<Real>(mesh->nverts() * ncomps * (dim + 1));
  auto f = LAMBDA(LO v) {
    if (!owned[v] || !interior[v]) return;
    auto qr_decomp = get_cavity_qr_decomposition<dim>(v, v2ve,
        ve2e, ev2v, coords);
    for (Int comp = 0; comp < ncomps; ++comp) {
      auto coeffs = fit_cavity_polynomial<dim>(qr_decomp, v,
          v2ve, ve2e, e_data, comp, ncomps);
      set_vector(out, v * ncomps + comp, coeffs);
    }
  };
  parallel_for(mesh->nverts(), f);
  return mesh->sync_array(VERT, Reals(out), ncomps * (dim + 1));
}

static Reals get_interior_coeffs(Mesh* mesh, Reals e_data, Int ncomps,
    Read<I8> interior) {
  if (mesh->dim() == 3) {
    return get_interior_coeffs_dim<3>(mesh, e_data, ncomps, interior);
  }
  if (mesh->dim() == 2) {
    return get_interior_coeffs_dim<2>(mesh, e_data, ncomps, interior);
  }
  NORETURN(Reals());
}

static void diffuse_to_exterior(Mesh* mesh, Reals* p_v_data, Int ncomps, Read<I8>* p_visited) {
  auto v_data = *p_v_data;
  auto visited = *p_visited;
  auto new_data = deep_copy(v_data);
  auto new_visited = deep_copy(visited);
  auto v2v = mesh->ask_star(VERT);
  auto v2vv = v2v.a2ab;
  auto vv2v = v2v.ab2b;
  auto f = LAMBDA(LO v) {
    if (visited[v]) return;
    Int nadj = 0;
    for (auto vv = v2vv[v]; vv < v2vv[v + 1]; ++vv) {
      auto ov = vv2v[vv];
      if (visited[ov]) ++nadj;
    }
    if (!nadj) return;
    for (Int comp = 0; comp < ncomps; ++comp) {
      Real sum = 0;
      for (auto vv = v2vv[v]; vv < v2vv[v + 1]; ++vv) {
        auto ov = vv2v[vv];
        if (visited[ov]) sum += v_data[ov * ncomps + comp];
      }
      new_data[v * ncomps + comp] = sum / nadj;
    }
    new_visited[v] = 1;
  };
  parallel_for(mesh->nverts(), f);
  v_data = new_data;
  visited = new_visited;
  v_data = mesh->sync_array(VERT, v_data, ncomps);
  visited = mesh->sync_array(VERT, visited, 1);
  *p_v_data = v_data;
  *p_visited = visited;
}

template <Int dim>
static Reals evaluate_coeffs_dim(Mesh* mesh, Reals v_coeffs, Int ncomps) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * ncomps);
  auto f = LAMBDA(LO v) {
    auto x = get_vector<dim>(coords, v);
    for (Int comp = 0; comp < ncomps; ++comp) {
      auto coeffs = get_vector<dim + 1>(v_coeffs, v * ncomps + comp);
      auto val = eval_polynomial(coeffs, x);
      out[v * ncomps + comp] = val;
    }
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals evaluate_coeffs(Mesh* mesh, Reals v_coeffs, Int ncomps) {
  if (mesh->dim() == 3) {
    return evaluate_coeffs_dim<3>(mesh, v_coeffs, ncomps);
  }
  if (mesh->dim() == 2) {
    return evaluate_coeffs_dim<2>(mesh, v_coeffs, ncomps);
  }
  NORETURN(Reals());
}

Reals project(Mesh* mesh, Reals e_data) {
  CHECK(mesh->owners_have_all_upward(VERT));
  CHECK(e_data.size() % mesh->nelems() == 0);
  auto ncomps = e_data.size() / mesh->nelems();
  auto class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto dim = mesh->dim();
  auto interior = each_eq_to(class_dim, I8(dim));
  auto have_local_interior = (max(interior) == 1);
  auto comm = mesh->comm();
  auto have_interior_verts = comm->reduce_or(have_local_interior);
  CHECK(have_interior_verts);
  auto v_coeffs = get_interior_coeffs(mesh, e_data, ncomps, interior);
  auto visited = interior;
  while (comm->reduce_or(min(visited) == 0)) {
    diffuse_to_exterior(mesh, &v_coeffs, ncomps * (dim + 1), &visited);
  }
  return evaluate_coeffs(mesh, v_coeffs, ncomps);
}

} // end namespace osh
