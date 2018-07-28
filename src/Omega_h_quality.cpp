#include "Omega_h_quality.hpp"

#include "Omega_h_for.hpp"

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
Reals measure_qualities_tmpl(Mesh* mesh, LOs a2e, Reals metrics) {
  MetricElementQualities<mesh_dim, metric_dim> measurer(mesh, metrics);
  auto ev2v = mesh->ask_verts_of(mesh_dim);
  auto na = a2e.size();
  Write<Real> qualities(na);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<mesh_dim + 1>(ev2v, e);
    qualities[a] = measurer.measure(v);
  };
  parallel_for(na, f, "measure_qualities");
  return qualities;
}

Reals measure_qualities(Mesh* mesh, LOs a2e, Reals metrics) {
  if (a2e.size() == 0) return Reals({});
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return measure_qualities_tmpl<3, 3>(mesh, a2e, metrics);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return measure_qualities_tmpl<2, 2>(mesh, a2e, metrics);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return measure_qualities_tmpl<3, 1>(mesh, a2e, metrics);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return measure_qualities_tmpl<2, 1>(mesh, a2e, metrics);
  }
  if (mesh->dim() == 1) {
    return Reals(a2e.size(), 1.0);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals measure_qualities(Mesh* mesh, LOs a2e) {
  return measure_qualities(mesh, a2e, mesh->get_array<Real>(VERT, "metric"));
}

Reals measure_qualities(Mesh* mesh) {
  return measure_qualities(mesh, LOs(mesh->nelems(), 0, 1));
}

Reals measure_qualities(Mesh* mesh, Reals metrics) {
  return measure_qualities(mesh, LOs(mesh->nelems(), 0, 1), metrics);
}

/* Cavity qualities are used for several things:
   1) Discarding operations that don't meet minimum requirements.
      In 1D, we will set our desired and allowed qualities to zero
      to ensure no operation is discarded for quality reasons
   2) Choosing among overlapping operations based on higher quality.
      This is the real key to Omega_h's performance.
      In 1D, the property of selecting highest quality is not so
      critical, because all qualities are essentially perfect.
      However, the special property that qualities must not form
      long monotonic paths in the conflict graph is at the center
      of Omega_h's ability to perform efficiently, especially in parallel.
      So we need to choose floating-point numbers such that no long
      monotonic chains are formed.
      This isn't too hard, what we'll essentially do is set values based
      on global numbers of keys being even or odd.
      This is of course assuming that the global numbering system keeps
      these numbers consecutive across X (which it should do).
      It must be based on global numbers to be deterministic in every
      sense that Omega_h is deterministic.
 */

Reals get_1d_cavity_qualities(Mesh* mesh, Int key_dim, LOs keys2kds) {
  auto kd_globals = mesh->globals(key_dim);
  auto nkeys = keys2kds.size();
  Write<Real> out(nkeys);
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto kd = keys2kds[key];
    auto global = kd_globals[kd];
    auto val = (global % 2);
    out[key] = Real(val);
  };
  parallel_for(nkeys, f, "get_1d_cavity_qualities");
  return out;
}

}  // end namespace Omega_h
