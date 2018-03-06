#include <Omega_h_solve.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_align.hpp>
#include <Omega_h_loop.hpp>

namespace Omega_h {

Reals matrix_vector_product(Mesh* mesh, Reals a_edge, Reals a_vert, Reals x) {
  auto verts2edges = mesh->ask_up(VERT, EDGE);
  auto edges2verts = mesh->ask_verts_of(EDGE);
  auto out_w = Write<Real>(x.size());
  auto f = OMEGA_H_LAMBDA(LO vert) {
    auto value = a_vert[vert] * x[vert];
    auto begin = verts2edges.a2ab[vert];
    auto end = verts2edges.a2ab[vert + 1];
    for (auto vert_edge = begin; vert_edge < end; ++vert_edge) {
      auto edge = verts2edges.ab2b[vert_edge];
      auto code = verts2edges.codes[vert_edge];
      auto edge_vert = code_which_down(code);
      auto other_vert = edges2verts[edge * 2 + (1 - edge_vert)];
      value += a_edge[edge] * x[other_vert];
    }
    out_w[vert] = value;
  };
  parallel_for(mesh->nverts(), f);
  return mesh->sync_array(VERT, Reals(out_w), 1);
}

Real vector_dot_product(CommPtr comm, Reals u, Reals v) {
  OMEGA_H_CHECK(u.size() == v.size());
  auto squares = multiply_each(u, v);
  return std::sqrt(repro_sum(comm, squares));
}

Reals conjugate_gradient(Mesh* mesh, Reals b,
    Reals a_edge, Reals a_vert, Reals x_0,
    Real tolerance, Int max_iters) {
  OMEGA_H_CHECK(tolerance >= 0.);
  auto comm = mesh->comm();
  auto x_k = x_0;
  auto ax_k = matrix_vector_product(mesh, a_edge, a_vert, x_0);
  auto r_k = subtract_each(b, ax_k);
  auto normsq_k = vector_dot_product(comm, r_k, r_k);
  if (normsq_k <= square(tolerance)) return x_k;
  auto p_k = r_k;
  for (Int k = 0; k < max_iters; ++k) {
    auto ap_k = matrix_vector_product(mesh, a_edge, a_vert, p_k);
    auto c_k = normsq_k / vector_dot_product(comm, p_k, ap_k);
    auto x_kp1 = add_each(x_k, multiply_each_by(p_k, c_k));
    auto r_kp1 = subtract_each(r_k, multiply_each_by(ap_k, c_k));
    auto normsq_kp1 = vector_dot_product(comm, r_kp1, r_kp1);
    if (normsq_kp1 <= square(tolerance)) return x_kp1;
    auto d_k = normsq_kp1 / normsq_k;
    auto p_kp1 = add_each(r_kp1, multiply_each_by(p_k, d_k));
    p_k = p_kp1;
    x_k = x_kp1;
    r_k = r_kp1;
    normsq_k = normsq_kp1;
  }
  Omega_h_fail("conjugate_gradient method failed to converge in %d iterations, residual norm is %.17e > %.17e\n",
      max_iters, std::sqrt(normsq_k), tolerance);
}

}
