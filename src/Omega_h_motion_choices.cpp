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
};

template <template <Int> typename Helper, Int dim>
MotionChoices motion_choices_tmpl(Mesh* mesh, AdaptOpts const& opts,
    LOs cands2verts) {
  using Metric = Helper<dim>::Metric;
  auto coords = mesh->coords();
  auto metrics = Helper<dim>::get_metrics_array(mesh);
  auto ncands = cands2verts.size();
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto max_steps = opts.max_motion_steps;
  CHECK(max_steps >= 0);
  auto step_size = opts.motion_step_size;
  auto max_length = opts.max_length_allowed;
  CHECK(0.0 < step_size);
  CHECK(step_size < 1.0);
  auto f = LAMBDA(LO cand) {
    auto v = cands2verts[cand];
    auto cm = Helper<dim>::get_metric(metrics, v);
    auto lcm = linearize_metric(mc);
    auto cx = get_vector<dim>(coords, v);
    bool first = true;
    Real last_qual = -1;
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
    }
  };
  parallel_for(ncands, f);
}

}
