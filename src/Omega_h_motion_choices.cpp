#include "Omega_h_motion.hpp"
#include "access.hpp"
#include "Omega_h_array_ops.hpp"
#include "loop.hpp"
#include "metric.hpp"
#include "quality.hpp"
#include "size.hpp"

namespace Omega_h {

template <Int dim>
class MetricMotion {
 public:
  using Metric = Matrix<dim, dim>;
  static constexpr Int ncomps = symm_dofs(dim);
  static Reals get_metrics_array(Mesh* mesh) {
    return mesh->get_array<Real>(VERT, "metric");
  }
  DEVICE static Metric get_metric(Reals metrics, LO v) {
    return get_symm<dim>(metrics, v);
  }
  template <typename Arr>
  DEVICE static Metric get_linear_metric(Arr const& tmp, Int offset) {
    Vector<symm_dofs(dim)> vec;
    for (Int i = 0; i < symm_dofs(dim); ++i) {
      vec[i] = tmp[offset + i];
    }
    return vector2symm(vec);
  }
  DEVICE static Metric gather_maxdet_metric(
      Reals metrics, Few<LO, dim + 1> kvv2v, Int kvv_c, Metric nm) {
    auto kvv2m = gather_symms<dim + 1, dim>(metrics, kvv2v);
    kvv2m[kvv_c] = nm;
    return maxdet_metric(kvv2m);
  }
};

template <Int dim>
class IsoMotion {
 public:
  using Metric = Real;
  static constexpr Int ncomps = 1;
  static Reals get_metrics_array(Mesh* mesh) {
    return mesh->get_array<Real>(VERT, "size");
  }
  DEVICE static Metric get_metric(Reals metrics, LO v) { return metrics[v]; }
  template <typename Arr>
  DEVICE static Metric get_linear_metric(Arr const& tmp, Int offset) {
    return tmp[offset];
  }
  /* metric is unused by shape quality measure, so don't bother
     gathering the other ones and actually computing some maximum */
  DEVICE static Metric gather_maxdet_metric(
      Reals, Few<LO, dim + 1>, Int, Metric nm) {
    return nm;
  }
};

template <Int dim>
INLINE Real metric_element_quality(Few<Vector<dim>, dim + 1> p, Real) {
  return metric_element_quality(p, DummyIsoMetric());
}

template <Int dim, typename Arr>
static DEVICE Vector<dim> get_coords(Arr const& tmp, Int offset) {
  Vector<dim> out;
  for (Int i = 0; i < dim; ++i) out[i] = tmp[offset + i];
  return out;
}

template <template <Int> class Helper, Int dim>
MotionChoices motion_choices_tmpl(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2verts) {
  using Metric = typename Helper<dim>::Metric;
  auto pack = pack_linearized_fields(mesh);
  constexpr Int maxcomps = 30;
  CHECK(pack.ncomps <= maxcomps);
  auto new_sol_w = deep_copy(pack.data);
  auto coords = mesh->coords();
  auto metrics = Helper<dim>::get_metrics_array(mesh);
  auto ncands = cands2verts.size();
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto v2k = mesh->ask_up(VERT, dim);
  auto kv2v = mesh->ask_verts_of(dim);
  auto verts2dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto cands2elems = unmap_graph(cands2verts, v2k);
  auto elems2old_qual = mesh->ask_qualities();
  auto cands2old_qual =
      graph_reduce(cands2elems, elems2old_qual, 1, OMEGA_H_MIN);
  auto max_steps = opts.max_motion_steps;
  CHECK(max_steps >= 0);
  auto step_size = opts.motion_step_size;
  auto max_length = opts.max_length_allowed;
  CHECK(0.0 < step_size);
  CHECK(step_size < 1.0);
  auto did_move_w = Write<I8>(ncands);
  auto coordinates_w = Write<Real>(ncands * dim);
  auto metrics_w = Write<Real>(ncands * Helper<dim>::ncomps);
  auto qualities_w = Write<Real>(ncands);
  auto f = LAMBDA(LO cand) {
    auto v = cands2verts[cand];
    auto v_dim = verts2dim[v];
    auto old_qual = cands2old_qual[cand];
    auto last_qual = old_qual;
    Vector<maxcomps> tmp;
    for (Int i = 0; i < pack.ncomps; ++i) {
      tmp[i] = pack.data[v * pack.ncomps + i];
    }
    for (Int step = 0; step < max_steps; ++step) {
      bool found_step = false;
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
          tmp[i] = (1.0 - step_size) * new_sol_w[v * pack.ncomps + i] +
                   step_size * pack.data[ov * pack.ncomps + i];
        }
        auto nx = get_coords<dim>(tmp, pack.coords_offset);
        auto lnm = Helper<dim>::get_linear_metric(tmp, pack.metric_offset);
        auto nm = delinearize_metric(lnm);
        auto overshoots = false;
        for (auto ve_l = v2e.a2ab[v]; ve_l < v2e.a2ab[v + 1]; ++ve_l) {
          auto e_l = v2e.ab2b[ve_l];
          auto ve_code_l = v2e.codes[ve_l];
          auto evv_c_l = code_which_down(ve_code_l);
          auto evv_o_l = 1 - evv_c_l;
          auto ov_l = ev2v[e_l * 2 + evv_o_l];
          auto om = Helper<dim>::get_metric(metrics, ov_l);
          auto ox = get_vector<dim>(coords, ov_l);
          Few<Vector<dim>, 2> evv2nx;
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
        for (auto vk = v2k.a2ab[v]; vk < v2k.a2ab[v + 1]; ++vk) {
          auto k = v2k.ab2b[vk];
          auto vk_code = v2k.codes[vk];
          auto kvv_c = code_which_down(vk_code);
          auto kvv2v = gather_verts<dim + 1>(kv2v, k);
          auto kvv2nx = gather_vectors<dim + 1, dim>(coords, kvv2v);
          kvv2nx[kvv_c] = nx;
          auto km =
              Helper<dim>::gather_maxdet_metric(metrics, kvv2v, kvv_c, nm);
          auto k_qual = metric_element_quality(kvv2nx, km);
          new_qual = min2(new_qual, k_qual);
          if (new_qual <= last_qual) break;
        }  // end loop over elements for quality
        if (new_qual <= last_qual) continue;
        found_step = true;
        last_qual = new_qual;
        for (Int i = 0; i < pack.ncomps; ++i) {
          new_sol_w[v * pack.ncomps + i] = tmp[i];
        }
      }  // end loop over edges for possible steps
      if (!found_step) break;
    }  // end loop over motion steps
    auto did_move = last_qual > old_qual;
    qualities_w[cand] = last_qual;
    did_move_w[cand] = did_move;
    if (did_move) CHECK(last_qual >= 0.0);
  };
  parallel_for(ncands, f);
  auto qualities = Reals(qualities_w);
  auto new_sol = Reals(new_sol_w);
  return {Read<I8>(did_move_w), Reals(qualities_w), Reals(new_sol_w)};
}

MotionChoices get_motion_choices(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2verts) {
  if (mesh->dim() == 3) {
    if (mesh->has_tag(VERT, "metric")) {
      return motion_choices_tmpl<MetricMotion, 3>(mesh, opts, cands2verts);
    } else {
      CHECK(mesh->has_tag(VERT, "size"));
      return motion_choices_tmpl<IsoMotion, 3>(mesh, opts, cands2verts);
    }
  } else {
    CHECK(mesh->dim() == 2);
    if (mesh->has_tag(VERT, "metric")) {
      return motion_choices_tmpl<MetricMotion, 2>(mesh, opts, cands2verts);
    } else {
      CHECK(mesh->has_tag(VERT, "size"));
      return motion_choices_tmpl<IsoMotion, 2>(mesh, opts, cands2verts);
    }
  }
}

}  // end namespace Omega_h
