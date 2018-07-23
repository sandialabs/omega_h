#ifndef OMEGA_H_QUALITY_HPP
#define OMEGA_H_QUALITY_HPP

#include <Omega_h_shape.hpp>

namespace Omega_h {

template <Int dim>
OMEGA_H_INLINE Real mean_ratio(Real size, Real mean_squared_length) {
  return power<2, dim>(size / equilateral_simplex_size(dim)) /
         mean_squared_length;
}

/* note that we will always use a constant metric tensor over the whole
 * element to compute its quality, because that way we are computing
 * the quality of the element after a single linear transformation which
 * is guaranteed not to turn it inside out.
 * this is why edge lengths are recomputed using the metric interpolated
 * to the element centroid.
 * other authors will try to use the highly accurate metric length
 * formula together with the very approximate metric volume formula to
 * compute quality. that can, for example, compute qualities greater
 * than 1.0, and other strange results.
 */

template <Int dim, typename Metric>
OMEGA_H_INLINE Real metric_element_quality(
    Few<Vector<dim>, dim + 1> p, Metric metric) {
  auto b = simplex_basis<dim, dim>(p);
  auto rs = simplex_size_from_basis(b);
  auto s = metric_size<dim>(rs, metric);
  if (s < 0) return s;
  auto ev = element_edge_vectors(p, b);
  auto msl = mean_squared_metric_length(ev, metric);
  return mean_ratio<dim>(s, msl);
}

template <Int space_dim, Int metric_dim>
struct MetricElementQualities {
  Reals coords;
  Reals metrics;
  MetricElementQualities(Mesh const* mesh, Reals metrics_in)
      : coords(mesh->coords()), metrics(metrics_in) {}
  MetricElementQualities(Mesh const* mesh)
      : MetricElementQualities(mesh, mesh->get_array<Real>(VERT, "metric")) {}
  OMEGA_H_DEVICE Real measure(Few<LO, space_dim + 1> v) const {
    auto p = gather_vectors<space_dim + 1, space_dim>(coords, v);
    auto ms = gather_symms<space_dim + 1, metric_dim>(metrics, v);
    auto m = maxdet_metric(ms);
    return metric_element_quality(p, m);
  }
};

Reals measure_qualities(Mesh* mesh, LOs a2e, Reals metrics);
Reals measure_qualities(Mesh* mesh, LOs a2e);
Reals measure_qualities(Mesh* mesh);
Reals measure_qualities(Mesh* mesh, Reals metrics);
Reals get_1d_cavity_qualities(Mesh* mesh, Int key_dim, LOs keys2kds);

}  // end namespace Omega_h

#endif
