#include "Omega_h_array_ops.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_motion.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_shape.hpp"

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
MotionChoices motion_choices_tmpl(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2verts) {
  auto tolerance = opts.motion_tolerance;
  OMEGA_H_CHECK(tolerance >= 0.0);
  auto nsteps = opts.motion_step_count;
  OMEGA_H_CHECK(nsteps > 0);
  auto step_size = opts.motion_step_size;
  OMEGA_H_CHECK(0.0 < step_size);
  OMEGA_H_CHECK(step_size < 1.0);
  auto coords = mesh->coords();
  auto metrics = mesh->get_array<Real>(VERT, "metric");
  auto size_errors = mesh->get_array<Real>(mesh_dim, "size_error");
  auto ncands = cands2verts.size();
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto v2k = mesh->ask_up(VERT, mesh_dim);
  auto kv2v = mesh->ask_verts_of(mesh_dim);
  auto verts2dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto cands2elems = unmap_graph(cands2verts, v2k);
  auto elems2old_qual = mesh->ask_qualities();
  auto cands2old_qual =
      graph_reduce(cands2elems, elems2old_qual, 1, OMEGA_H_MIN);
  auto new_coords_w = Write<Real>(ncands * mesh_dim);
  auto new_quals_w = Write<Real>(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto v = cands2verts[cand];
    Real old_obj = -1.0;
    auto old_x = get_vector<mesh_dim>(coords, v);
    auto new_x = old_x;
    bool did_move = false;
    while (true) {
      Real new_obj = 0.0;
      auto obj_grad = zero_vector<mesh_dim>();
      bool obj_converged = true;
      bool quality_ok = true;
      for (auto vk = v2k.a2ab[v]; vk < v2k.a2ab[v + 1]; ++vk) {
        auto k = v2k.ab2b[vk];
        auto vk_code = v2k.codes[vk];
        auto kvv_c = code_which_down(vk_code);
        auto kvv2v = gather_verts<mesh_dim + 1>(kv2v, k);
        auto kvv2nx = gather_vectors<mesh_dim + 1, mesh_dim>(coords, kvv2v);
        kvv2nx[kvv_c] = new_x;
        auto k_basis = simplex_basis(kvv2nx);
        auto k_tmp_size = element_size(k_basis);
        auto k_size_grad = get_size_gradient(kvv2nx, kvv_c);
        auto k_size_diff = k_tmp_size - orig_sizes[k];
        auto k_tmp_error = size_errors[k] + k_size_diff;
        auto k_tmp_rel_error = k_tmp_error / k_tmp_size;
        if (k_tmp_rel_error > tolerance) obj_converged = false;
        auto k_obj = square(k_tmp_error);
        auto k_obj_grad = 2.0 * k_tmp_error * k_size_grad;
        auto kvv2m = gather_symms<mesh_dim + 1, metric_dim>(metrics, kvv2v);
        auto km = maxdet_metric(kvv2m);
        auto k_qual = metric_element_quality(kvv2nx, km);
        if (k_qual < min_qual_allowed) {
          quality_ok = false;
          break;
        }
        new_obj += k_obj;
        obj_grad += k_obj_grad;
      }
      if (step && (!quality_ok || new_obj >= old_obj))
        /* we either took a step or backtracked, and the new position
           is bad either because it increases the objective or the elements
           are too distorted */
        if (step == nsteps) {
          /* out of time, fully back up to the last good position */
          new_x = old_x;
          break;
        } else {
          /* still have time, backtrack */
          new_x = (new_x + old_x) / 2.0;
        }
      } else {
        /* either no prior step exists, or the prior step both reduced
           the objective (good) and created acceptable quality elements */
        if (step) did_move = true;
        if (obj_converged || step == nsteps) {
          /* either we've solved the problem or we're out of time.
             the last step (if any) was a good one. don't move. */
          break;
        } else {
          /* actually take a new step */
          old_x = new_x;
          old_obj = new_obj;
          new_x = (old_x - (obj_grad * step_size));
        }
      }
    }
    did_move_w[cand] = I8(did_move);
    new_quals_w[cand] = new_qual;
    set_vector(new_coords_w, cand, new_x);
  };
  parallel_for(ncands, f);
  auto new_quals = Reals(new_quals_w);
  auto new_coords = Reals(new_coords_w);
  return {new_quals, new_coords};
}

MotionChoices get_motion_choices(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2verts) {
  auto metric_dim = get_metric_dim(mesh);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return motion_choices_tmpl<3, 3>(mesh, opts, cands2verts);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return motion_choices_tmpl<2, 2>(mesh, opts, cands2verts);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return motion_choices_tmpl<3, 1>(mesh, opts, cands2verts);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return motion_choices_tmpl<2, 1>(mesh, opts, cands2verts);
  }
  OMEGA_H_NORETURN(MotionChoices());
}

}  // end namespace Omega_h
