#include "Omega_h_motion.hpp"
#include "metric.hpp"
#include "access.hpp"

namespace Omega_h {

template <Int dim>
class MetricMotion<dim> {
  using Metric = Matrix<dim, dim>;
  static Reals get_metrics_array(Mesh* mesh) {
    return mesh->get_array<Real>(VERT, "metric");
  }
  DEVICE static Metric get_metric(Reals metrics, LO v) {
    return get_symm<dim>(metrics, v);
  }
  DEVICE static Metric gather_maxdet_metric(Few<LO, dim + 1> kvv2v,
      Int kvv_c, Metric nm) {
    auto kvv2m = gather_symms<dim + 1, dim>(metrics, kvv2v);
    kvv2m[kvv_c] = nm;
    return maxdet_metric(kvv2m);
  }
};

template <Int dim>
class IsoMotion<dim> {
  using Metric = Real;
  static Reals get_metrics_array(Mesh* mesh) {
    return mesh->get_array<Real>(VERT, "size");
  }
  DEVICE static Metric get_metric(Reals metrics, LO v) {
    return metrics[v];
  }
  /* metric is unused by shape quality measure, so don't bother
     gathering the other ones and actually computing some maximum */
  DEVICE static Metric gather_maxdet_metric(Few<LO, dim + 1>,
      Int, Metric nm) {
    return nm;
  }
};

template <Int dim>
INLINE Real metric_element_quality(Few<Vector<dim>, dim + 1> p, Real h) {
  return metric_element_quality(p, DummyIsoMetric());
}

template <template <Int> typename Helper, Int dim>
MotionChoices motion_choices_tmpl(Mesh* mesh, AdaptOpts const& opts,
    LOs cands2verts) {
  using Metric = Helper<dim>::Metric;
  auto coords = mesh->coords();
  auto metrics = Helper<dim>::get_metrics_array(mesh);
  auto ncands = cands2verts.size();
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto v2k = mesh->ask_up(VERT, dim);
  auto cands2elems = unmap_graph(cands2verts, v2k);
  auto elems2old_qual = mesh->ask_qualities();
  auto cands2old_qual = graph_reduce(cands2elems, elems2old_qual, 1, OMEGA_H_MIN);
  auto max_steps = opts.max_motion_steps;
  CHECK(max_steps >= 0);
  auto step_size = opts.motion_step_size;
  auto max_length = opts.max_length_allowed;
  CHECK(0.0 < step_size);
  CHECK(step_size < 1.0);
  auto coordinates_w = Write<Real>(ncands * dim);
  auto qualities_w = Write<Real>(ncands);
  auto f = LAMBDA(LO cand) {
    auto v = cands2verts[cand];
    auto cm = Helper<dim>::get_metric(metrics, v);
    auto lcm = linearize_metric(mc);
    auto cx = get_vector<dim>(coords, v);
    Real last_qual = cands2old_qual[v];
    for (Int step = 0; step < max_steps; ++step) {
      auto best_qual = last_qual;
      bool found_step = false;
      auto best_x = cx;
      auto best_m = cm;
      for (auto ve = v2e.a2ab[v]; ve < v2e.a2ab[v + 1]; ++ve) {
        auto e = v2e.ab2b[ve];
        auto evv_c = v2e.codes[ve];
        auto evv_o = 1 - evv_c;
        auto evv2v = gather_verts<2>(ev2v, e);
        auto ov = evv2v[evv_o];
        auto ox = get_vector<dim>(coords, ov);
        auto d = ox - cx;
        auto s = d * step_size;
        auto nx = cx + s;
        auto om = Helper<dim>::get_metric(metrics, ov);
        auto lom = linearize_metric(om); // could cache these
        auto lnm = lcm * (1.0 - step_size) + lom * step_size;
        auto nm = delinearize_metric(lnm);
        Few<Vector<dim> 2> evv2nx;
        Few<Metric, 2> evv2nm;
        evv2nx[evv_c] = nx;
        evv2nx[evv_o] = ox;
        evv2nm[evv_c] = nm;
        evv2nm[evv_o] = om;
        auto nl = metric_edge_length(evv2nx, evv2nm);
        if (nl > max_length) continue;
        Real new_qual = 1.0;
        for (auto vk = v2k.a2ab[v]; vk < v2k.a2ab[v + 1]; ++vk) {
          auto k = v2k.ab2b[vk];
          auto kvv_c = v2k.codes[vk];
          auto kvv2v = gather_verts<dim + 1>(k2v, k);
          auto kvv2nx = gather_vectors<dim + 1, dim>(coords, kvv2v);
          kvv2nx[kvv_c] = nx;
          auto km = Helper<dim>::gather_maxdet_metric(kvv2v, kvv_c, nm);
          auto k_qual = metric_element_quality(kvv2nx, km);
          new_qual = min2(new_qual, k_qual);
          if (new_qual <= best_qual) break;
        } // end loop over elements for quality
        if (new_qual <= best_qual) continue;
        found_step = true;
        best_qual = new_qual;
        best_x = nx;
        best_m = nm;
      } // end loop over edges for possible steps
      if (!found_step) break;
      last_qual = best_qual;
      cx = best_x;
      cm = best_m;
      lcm = linearize_metric(mc);
    } // end loop over motion steps
    set_vector(coordinates_w, cand, cx);
    qualities_w[cand] = cm;
  };
  parallel_for(ncands, f);
  return { Reals(qualities_w), Reals(coordinates_w) };
}

MotionChoices get_motion_choices(Mesh* mesh, AdaptOpts const& opts,
    LOs cands2verts) {
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

}
