#include "Omega_h_motion.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_shape.hpp"

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
MotionChoices motion_choices_tmpl(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2verts) {
  auto max_steps = opts.xfer_opts.max_size_steps;
  OMEGA_H_CHECK(max_steps > 0);
  auto min_step_size = opts.xfer_opts.min_size_step_ratio;
  OMEGA_H_CHECK(0.0 <= min_step_size);
  OMEGA_H_CHECK(min_step_size <= 1.0);
  auto max_rel_error = opts.xfer_opts.max_size_error_ratio;
  OMEGA_H_CHECK(max_rel_error >= 0.0);
  auto min_qual_allowed = opts.min_quality_allowed;
  auto max_ml_allowed = opts.max_length_allowed;
  auto coords = mesh->coords();
  auto metrics = mesh->get_array<Real>(VERT, "metric");
  auto size_errors = mesh->get_array<Real>(mesh_dim, "size_error");
  auto orig_sizes = mesh->ask_sizes();
  auto ncands = cands2verts.size();
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto v2k = mesh->ask_up(VERT, mesh_dim);
  auto kv2v = mesh->ask_verts_of(mesh_dim);
  auto did_move_w = Write<I8>(ncands);
  auto new_coords_w = Write<Real>(ncands * mesh_dim);
  auto new_quals_w = Write<Real>(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto v = cands2verts[cand];
    Real old_obj = -1.0;
    auto orig_x = get_vector<mesh_dim>(coords, v);
    auto old_x = orig_x;
    auto new_x = old_x;
    bool did_move = false;
    Real new_qual;
    Real max_rl = 0.0;
    for (auto ve = v2e.a2ab[v]; ve < v2e.a2ab[v + 1]; ++ve) {
      auto e = v2e.ab2b[ve];
      auto ve_code = v2e.codes[ve];
      auto eev_c = code_which_down(ve_code);
      auto ov = ev2v[e * 2 + (1 - eev_c)];
      auto ox = get_vector<mesh_dim>(coords, ov);
      auto rl = norm(orig_x - ox);
      max_rl = max2(max_rl, rl);
    }
    for (Int step = 0; true; ++step) {
      Real new_obj = 0.0;
      auto obj_grad = zero_vector<mesh_dim>();
      bool obj_converged = true;
      bool quality_ok = true;
      new_qual = 1.0;
      for (auto vk = v2k.a2ab[v]; vk < v2k.a2ab[v + 1]; ++vk) {
        auto k = v2k.ab2b[vk];
        auto vk_code = v2k.codes[vk];
        auto kkv_c = code_which_down(vk_code);
        auto kkv2v = gather_verts<mesh_dim + 1>(kv2v, k);
        auto kkv2nx = gather_vectors<mesh_dim + 1, mesh_dim>(coords, kkv2v);
        kkv2nx[kkv_c] = new_x;
        auto k_basis = simplex_basis<mesh_dim, mesh_dim>(kkv2nx);
        auto k_tmp_size = simplex_size_from_basis(k_basis);
        auto k_size_grad = get_size_gradient(kkv2nx, kkv_c);
        auto k_size_diff = k_tmp_size - orig_sizes[k];
        auto k_tmp_error = size_errors[k] + k_size_diff;
        auto k_tmp_rel_error = std::abs(k_tmp_error / k_tmp_size);
        if (k_tmp_rel_error > max_rel_error) {
          obj_converged = false;
        }
        auto k_obj = square(k_tmp_error);
        auto k_obj_grad = 2.0 * k_tmp_error * k_size_grad;
        auto kvv2m = gather_symms<mesh_dim + 1, metric_dim>(metrics, kkv2v);
        auto km = maxdet_metric(kvv2m);
        auto k_qual = metric_element_quality(kkv2nx, km);
        if (k_qual < min_qual_allowed) {
          quality_ok = false;
          break;
        }
        new_obj += k_obj;
        obj_grad += k_obj_grad;
        new_qual = min2(new_qual, k_qual);
      }
      Real max_ml = 0.0;
      for (auto ve = v2e.a2ab[v]; ve < v2e.a2ab[v + 1]; ++ve) {
        auto e = v2e.ab2b[ve];
        auto ve_code = v2e.codes[ve];
        auto eev_c = code_which_down(ve_code);
        auto eev2v = gather_verts<2>(ev2v, e);
        auto eev2nx = gather_vectors<2, mesh_dim>(coords, eev2v);
        eev2nx[eev_c] = new_x;
        auto eev2m = gather_symms<2, metric_dim>(metrics, eev2v);
        auto ml = metric_edge_length(eev2nx, eev2m);
        max_ml = max2(max_ml, ml);
      }
      if (step &&
          (!quality_ok || new_obj >= old_obj || max_ml > max_ml_allowed)) {
        /* we either took a step or backtracked, and the new position
           is bad either because it increases the objective, or the elements
           are too distorted, or an edge got too long */
        if (step == max_steps) {
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
        if (obj_converged || step == max_steps) {
          /* either we've solved the problem or we're out of time.
             the last step (if any) was a good one. don't move. */
          break;
        } else {
          /* actually take a new step */
          old_x = new_x;
          old_obj = new_obj;
          Vector<mesh_dim> delta_x;
          for (Int i = 0; i < mesh_dim; ++i) {
            delta_x[i] = -(new_obj / obj_grad[i]);
          }
          new_x = old_x + delta_x;
        }
      }
    }
    if (did_move && (norm(new_x - orig_x) < (min_step_size * max_rl))) {
      new_x = orig_x;
      did_move = false;
    }
    did_move_w[cand] = I8(did_move);
    new_quals_w[cand] = new_qual;
    set_vector(new_coords_w, cand, new_x);
  };
  parallel_for(ncands, f);
  auto did_move = Bytes(did_move_w);
  auto new_quals = Reals(new_quals_w);
  auto new_coords = Reals(new_coords_w);
  did_move = mesh->sync_subset_array(VERT, did_move, cands2verts, I8(-1), 1);
  new_quals = mesh->sync_subset_array(VERT, new_quals, cands2verts, -1.0, 1);
  new_coords =
      mesh->sync_subset_array(VERT, new_coords, cands2verts, 0.0, mesh_dim);
  return {did_move, new_quals, new_coords};
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
