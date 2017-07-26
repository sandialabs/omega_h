#include "Omega_h_array_ops.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_motion.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_shape.hpp"

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
struct MotionChooser {
 private:
  static constexpr Int metric_ncomps = symm_ncomps(metric_dim);
 private:
  LinearPack pack;
  Graph v2k;
  LOs kv2v;
  Reals coords;
  Reals metrics;
  /* this array contains the element size (volume) errors as they are
     at the beginning of each step */
  Write<Real> size_errors_w;
  /* this array contains the vertex field DOFs as they are
     at the beginning of each step */
  Write<Real> new_sol_w;
  /* this array contains the vertex field DOFs as they are
     at each candidate point to be evaluated during a step */
  Write<Real> tmp_sol_w;
 private:
  OMEGA_H_INLINE Real compute_objective(LO v) const {
    Real obj = 0.0;
    auto nx = get_interpolated_position(v);
    for (auto vk = v2k.a2ab[v]; vk < v2k.a2ab[v + 1]; ++vk) {
      auto k = v2k.ab2b[vk];
      auto k_obj = size_errors_w[k];
      obj += square(k_obj);
    }  // end loop over elements for quality
    return obj;
  }
  OMEGA_H_INLINE Vector<mesh_dim> compute_objective_gradient(LO v) const {
    auto obj_grad = zero_vector<mesh_dim>();
    auto nx = get_tmp_position(v);
    for (auto vk = v2k.a2ab[v]; vk < v2k.a2ab[v + 1]; ++vk) {
      auto k = v2k.ab2b[vk];
      auto vk_code = v2k.codes[vk];
      auto kvv_c = code_which_down(vk_code);
      auto kvv2v = gather_verts<mesh_dim + 1>(kv2v, k);
      auto kvv2nx = gather_vectors<mesh_dim + 1, mesh_dim>(coords, kvv2v);
      kvv2nx[kvv_c] = nx;
      auto k_size_grad = get_size_gradient(kvv2nx, kvv_c);
      auto k_obj = size_errors_w[k];
      auto k_obj_grad = 2.0 * k_obj * k_size_grad;
      obj_grad += k_obj_grad;
    }
    return obj_grad;
  }
  OMEGA_H_INLINE Real compute_quality(LO v) const {
    Real new_qual = 1.0;
    auto nx = get_tmp_position(v);
    auto nm = delinearize_metric(get_tmp_linear_metric(v));
    for (auto vk = v2k.a2ab[v]; vk < v2k.a2ab[v + 1]; ++vk) {
      auto k = v2k.ab2b[vk];
      auto vk_code = v2k.codes[vk];
      auto kvv_c = code_which_down(vk_code);
      auto kvv2v = gather_verts<mesh_dim + 1>(kv2v, k);
      auto kvv2nx = gather_vectors<mesh_dim + 1, mesh_dim>(coords, kvv2v);
      kvv2nx[kvv_c] = nx;
      auto kvv2m = gather_symms<mesh_dim + 1, metric_dim>(metrics, kvv2v);
      kvv2m[kvv_c] = nm;
      auto km = maxdet_metric(kvv2m);
      auto k_qual = metric_element_quality(kvv2nx, km);
      new_qual = min2(new_qual, k_qual);
    }  // end loop over elements for quality
    return new_qual;
  }
  OMEGA_H_INLINE Vector<mesh_dim> get_tmp_position(LO v) {
    Vector<mesh_dim> nx;
    for (Int i = 0; i < mesh_dim; ++i) {
      nx[i] = tmp_sol_w[v * pack.ncomps + pack.coords_offset + i];
    }
    return nx;
  }
  OMEGA_H_INLINE Matrix<metric_dim, metric_dim> get_tmp_linear_metric(LO v) {
    Vector<metric_ncomps> metric_comps;
    for (Int i = 0; i < metric_ncomps; ++i) {
      metric_comps[i] = tmp_sol_w[v * pack.ncomps + pack.metric_offset + i];
    }
    return vector2symm(metric_comps);
  }
};

template <Int mesh_dim, Int metric_dim>
MotionChoices motion_choices_tmpl(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2verts) {
  using Metric = Matrix<metric_dim, metric_dim>;
  auto pack = pack_linearized_fields(mesh, opts.xfer_opts);
  /* this array contains the vertex field DOFs as they are
     at the beginning of each step */
  auto new_sol_w = deep_copy(pack.data);
  /* this array contains the vertex field DOFs as they are
     at each candidate point to be evaluated during a step */
  auto tmp_sol_w = deep_copy(pack.data);
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
  auto max_steps = opts.max_motion_steps;
  OMEGA_H_CHECK(max_steps >= 0);
  auto max_backtracks = opts.max_motion_backtracks;
  auto max_length = opts.max_length_allowed;
  OMEGA_H_CHECK(0.0 < step_size);
  OMEGA_H_CHECK(step_size < 1.0);
  auto did_move_w = Write<I8>(ncands);
  auto coordinates_w = Write<Real>(ncands * mesh_dim);
  constexpr auto metric_ncomps = symm_ncomps(metric_dim);
  auto metrics_w = Write<Real>(ncands * metric_ncomps);
  auto qualities_w = Write<Real>(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto v = cands2verts[cand];
    auto v_dim = verts2dim[v];
    auto old_qual = cands2old_qual[cand];
    auto last_qual = old_qual;
    // loop over number of steps */
    for (Int step = 0; step < max_steps; ++step) {
      bool found_step = false;
      // loop over edges for possible steps
      for (auto ve = v2e.a2ab[v]; ve < v2e.a2ab[v + 1]; ++ve) {
        auto e = v2e.ab2b[ve];
        auto e_dim = edges2dim[e];
        if (e_dim != v_dim) continue;
        auto ve_code = v2e.codes[ve];
        auto evv_c = code_which_down(ve_code);
        auto evv_o = 1 - evv_c;
        auto evv2v = gather_verts<2>(ev2v, e);
        auto ov = evv2v[evv_o];
        for (Int i = 0; i < pack.ncomps; ++i) {
          tmp_sol_w[v * pack.ncomps + i] =
            (1.0 - step_size) * new_sol_w[v * pack.ncomps + i] +
                    step_size * pack.data[ov * pack.ncomps + i];
        }
        Vector<mesh_dim> nx;
        for (Int i = 0; i < mesh_dim; ++i) {
          nx[i] = tmp_sol_w[v * pack.ncomps + pack.coords_offset + i];
        }
        Vector<metric_ncomps> metric_comps;
        for (Int i = 0; i < metric_ncomps; ++i) {
          metric_comps[i] = tmp_sol_w[v * pack.ncomps + pack.metric_offset + i];
        }
        auto lnm = vector2symm(metric_comps);
        auto nm = delinearize_metric(lnm);
        auto overshoots = false;
        for (auto ve_l = v2e.a2ab[v]; ve_l < v2e.a2ab[v + 1]; ++ve_l) {
          auto e_l = v2e.ab2b[ve_l];
          auto ve_code_l = v2e.codes[ve_l];
          auto evv_c_l = code_which_down(ve_code_l);
          auto evv_o_l = 1 - evv_c_l;
          auto ov_l = ev2v[e_l * 2 + evv_o_l];
          auto om = get_symm<metric_dim>(metrics, ov_l);
          auto ox = get_vector<mesh_dim>(coords, ov_l);
          Few<Vector<mesh_dim>, 2> evv2nx;
          Few<Metric, 2> evv2nm;
          evv2nx[evv_c_l] = nx;
          evv2nx[evv_o_l] = ox;
          evv2nm[evv_c_l] = nm;
          evv2nm[evv_o_l] = om;
          auto nl = metric_edge_length(evv2nx, evv2nm);
          if (nl > max_length) {
            overshoots = true;
            break;
          }
        }  // end loop over edges again for overshooting
        if (overshoots) continue;
        Real new_qual = 1.0;
        // loop over elements to check new qualities
        for (auto vk = v2k.a2ab[v]; vk < v2k.a2ab[v + 1]; ++vk) {
          auto k = v2k.ab2b[vk];
          auto vk_code = v2k.codes[vk];
          auto kvv_c = code_which_down(vk_code);
          auto kvv2v = gather_verts<mesh_dim + 1>(kv2v, k);
          auto kvv2nx = gather_vectors<mesh_dim + 1, mesh_dim>(coords, kvv2v);
          kvv2nx[kvv_c] = nx;
          auto kvv2m = gather_symms<mesh_dim + 1, metric_dim>(metrics, kvv2v);
          kvv2m[kvv_c] = nm;
          auto km = maxdet_metric(kvv2m);
          auto k_qual = metric_element_quality(kvv2nx, km);
          new_qual = min2(new_qual, k_qual);
          if (new_qual <= last_qual) break;
        }  // end loop over elements for quality
        if (new_qual <= last_qual) continue;
        found_step = true;
        last_qual = new_qual;
        for (Int i = 0; i < pack.ncomps; ++i) {
          new_sol_w[v * pack.ncomps + i] = tmp_sol_w[v * pack.ncomps + i];
        }
      }  // end loop over edges for possible steps
      if (!found_step) break;
    }  // end loop over motion steps
    auto did_move = last_qual > old_qual;
    qualities_w[cand] = last_qual;
    did_move_w[cand] = did_move;
    if (did_move) OMEGA_H_CHECK(last_qual >= 0.0);
  };
  parallel_for(ncands, f);
  auto qualities = Reals(qualities_w);
  auto new_sol = Reals(new_sol_w);
  return {Read<I8>(did_move_w), Reals(qualities_w), Reals(new_sol_w)};
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
