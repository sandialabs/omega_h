#ifndef QUALITY_HPP
#define QUALITY_HPP

#include "algebra.hpp"
#include "host_few.hpp"
#include "size.hpp"

namespace Omega_h {

template <Int dim>
struct EquilateralSize;

template <>
struct EquilateralSize<2> {
  static constexpr Real value = 0.4330127018922193;  // sqrt(3)/4
};

template <>
struct EquilateralSize<3> {
  static constexpr Real value = 0.1178511301977579;  // 1/sqrt(72)
};

template <Int dim>
INLINE constexpr Real equilateral_size() {
  return EquilateralSize<dim>::value;
}

template <Int dim>
INLINE Real mean_ratio(Real size, Real mean_squared_length) {
  return power<2, dim>(size / equilateral_size<dim>()) / mean_squared_length;
}

INLINE Real metric_size(Real real_size, NoMetric) { return real_size; }

/* This paper (and a few others):
 *
 * Loseille, Adrien, Victorien Menier, and Frederic Alauzet.
 * "Parallel Generation of Large-size Adapted Meshes."
 * Procedia Engineering 124 (2015): 57-69.
 *
 * Mentions using $\sqrt{\det(M)}$ to compute volume in metric space.
 *
 * The call to power() allows us to pass in a 1x1 isotropic metric,
 * and have its "determinant" raised to the right power before taking
 * the square root, and even accepting its existing value in the case
 * of space being 2D.
 */

template <Int space_dim, Int metric_dim>
INLINE Real metric_size(Real real_size, Matrix<metric_dim, metric_dim> metric) {
  return real_size * power<space_dim, 2 * metric_dim>(determinant(metric));
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
INLINE Real metric_element_quality(Few<Vector<dim>, dim + 1> p, Metric metric) {
  auto b = simplex_basis<dim, dim>(p);
  auto rs = element_size(b);
  auto s = metric_size(rs, metric);
  if (s < 0) return s;
  auto ev = element_edge_vectors(p, b);
  auto msl = mean_squared_metric_length(ev, metric);
  return mean_ratio<dim>(s, msl);
}

template <Int dim>
INLINE Real real_element_quality(Few<Vector<dim>, dim + 1> p) {
  return metric_element_quality(p, NoMetric());
}

struct RealElementQualities {
  Reals coords;
  RealElementQualities(Mesh const* mesh) : coords(mesh->coords()) {}
  template <Int space_dim>
  DEVICE Real measure(Few<LO, space_dim + 1> v) const {
    auto p = gather_vectors<space_dim + 1, space_dim>(coords, v);
    return real_element_quality(p);
  }
};

struct MetricElementQualities {
  Reals coords;
  Reals metrics;
  MetricElementQualities(Mesh const* mesh)
      : coords(mesh->coords()),
        metrics(mesh->get_array<Real>(VERT, "metric")) {}
  template <Int space_dim, Int metric_dim>
  DEVICE Real measure(Few<LO, space_dim + 1> v) const {
    auto p = gather_vectors<space_dim + 1, space_dim>(coords, v);
    auto ms = gather_symms<space_dim + 1, metric_dim>(metrics, v);
    auto m = maxdet_metric(ms);
    return metric_element_quality(p, m);
  }
};

Reals measure_qualities(Mesh* mesh, LOs a2e);
Reals measure_qualities(Mesh* mesh);

}  // end namespace Omega_h

#endif
