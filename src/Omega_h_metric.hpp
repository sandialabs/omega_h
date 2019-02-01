#ifndef OMEGA_H_METRIC_HPP
#define OMEGA_H_METRIC_HPP

#include <Omega_h_eigen.hpp>

namespace Omega_h {

template <Int dim>
OMEGA_H_INLINE Real metric_product(Tensor<dim> const m, Vector<dim> const v) {
  return v * (m * v);
}

template <Int space_dim>
OMEGA_H_INLINE typename std::enable_if<(space_dim > 1), Real>::type
metric_product(Tensor<1> const m, Vector<space_dim> const v) {
  return v * (m[0][0] * v);
}

template <Int metric_dim, Int space_dim>
OMEGA_H_INLINE Real metric_length(
    Tensor<metric_dim> const m, Vector<space_dim> const v) {
  return std::sqrt(metric_product(m, v));
}

template <Int dim>
OMEGA_H_INLINE Real metric_desired_length(
    Tensor<dim> const m, Vector<dim> const dir) {
  return 1.0 / metric_length(m, dir);
}

OMEGA_H_INLINE Real metric_length_from_eigenvalue(Real const l) {
  return 1.0 / std::sqrt(l);
}

OMEGA_H_INLINE Real metric_eigenvalue_from_length(Real const h) {
  return 1.0 / square(h);
}

template <Int dim>
OMEGA_H_INLINE Vector<dim> metric_lengths_from_eigenvalues(
    Vector<dim> const l) {
  Vector<dim> h;
  for (Int i = 0; i < dim; ++i) h[i] = metric_length_from_eigenvalue(l[i]);
  return h;
}

template <Int dim>
OMEGA_H_INLINE Vector<dim> metric_eigenvalues_from_lengths(
    Vector<dim> const h) {
  Vector<dim> l;
  for (Int i = 0; i < dim; ++i) l[i] = metric_eigenvalue_from_length(h[i]);
  return l;
}

template <Int dim>
OMEGA_H_INLINE Tensor<dim> compose_metric(
    Tensor<dim> const r, Vector<dim> const h) {
  auto const l = metric_eigenvalues_from_lengths(h);
  return compose_ortho(r, l);
}

template <Int dim>
OMEGA_H_INLINE_BIG DiagDecomp<dim> decompose_metric(Tensor<dim> const m) {
  auto const ed = decompose_eigen(m);
  auto const h = metric_lengths_from_eigenvalues(ed.l);
  return {ed.q, h};
}

/* Alauzet details four different ways to interpolate
   the metric tensor:

1) M(t) = ((1-t)M_1^{-1/2} + t M_2^{-1/2})^{-2}

2) M(t) = (M_1^{-1/2} (M_2^{-1/2} / M_1^{-1/2})^t)^2

3) M(t) = ((1-t)M_1^{-1} + tM_2^{-1})^{-1}

4) M(t) = (1-t)M_1 + t M_2

The first three exhibit decent interpolation behavior.
The last one, per-component linear interpolation,
tends to produce very small isotropic ellipsoids given
two anisotropic ellipsoids, so is not good enough.
Both (1) and (2) require an eigendecomposition to get M_i^{-1/2},
which is relatively expensive.
Both (2) and (3) can be generalized to multiple input
tensors, for interpolation in a triangle or tet.

Looking a (1), (2) and (3) suggests that their only
difference is an operation we will call "linearization",
in this case converting the metric tensor into a quantity
that can be safely linearly interpolated.
(1) M^{-1/2}
(2) M^{-1}
(3) M

There is a fifth (fourth ?) option advocated by Loseille,
Michal, and Krakos which is to use the matrix logarithm
of M as the "linearized" quantity.
This is also consistent with work by Mota on using Lie
algebras to interpolate tensor quantities.
That is the mechanism we use here:
*/

template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> linearize_metric(Tensor<dim> const m) {
  return log_spd_old(m);
}

template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> delinearize_metric(Tensor<dim> const log_m) {
  return exp_spd_old(log_m);
}

template <Int n, typename T>
OMEGA_H_INLINE_BIG Few<T, n> linearize_metrics(Few<T, n> ms) {
  Few<T, n> log_ms;
  for (Int i = 0; i < n; ++i) log_ms[i] = linearize_metric(ms[i]);
  return log_ms;
}

/* the "proper" way to interpolate the metric tensor to
 * the barycenter of a simplex; does several eigendecompositions
 */
template <Int dim>
OMEGA_H_INLINE_BIG void average_metric_contrib(
    Tensor<dim>& am, Int& n, Tensor<dim> const m, bool const has_degen) {
  if (has_degen && max_norm(m) < OMEGA_H_EPSILON) return;
  am += linearize_metric(m);
  n++;
}

template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> average_metric_finish(
    Tensor<dim> am, Int const n, bool const has_degen) {
  if (has_degen && n == 0) return am;
  am /= n;
  return delinearize_metric(am);
}

template <Int dim, Int n>
OMEGA_H_INLINE_BIG Tensor<dim> average_metric(
    Few<Tensor<dim>, n> const ms, bool const has_degen) {
  auto am = zero_matrix<dim, dim>();
  Int ngood = 0;
  for (Int i = 0; i < n; ++i) {
    average_metric_contrib(am, ngood, ms[i], has_degen);
  }
  return average_metric_finish(am, ngood, has_degen);
}

template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> clamp_metric(
    Tensor<dim> const m, Real const h_min, Real const h_max) {
  auto ed = decompose_eigen(m);
  auto const l_max = metric_eigenvalue_from_length(h_min);
  auto const l_min = metric_eigenvalue_from_length(h_max);
  for (Int i = 0; i < dim; ++i) ed.l[i] = clamp(ed.l[i], l_min, l_max);
  return compose_ortho(ed.q, ed.l);
}

/* a cheap hackish variant of interpolation for getting a metric
 * tensor to use to measure an element's quality.
 * basically, choose the one that is asking for the smallest real-space volume
 * (big determinant means large metric volume which triggers refinement)
 * the reason we use a cheap hack is because the Log-Euclidean interpolation
 * we use is rather expensive, and we'd like to avoid calling it for every
 * potential element (we do a lot of cavity pre-evaluation).
 */
template <Int dim, Int n>
OMEGA_H_INLINE_BIG Tensor<dim> maxdet_metric(Few<Tensor<dim>, n> const ms) {
  auto m = ms[0];
  auto maxdet = determinant(m);
  for (Int i = 1; i < n; ++i) {
    auto det = determinant(ms[i]);
    if (det > maxdet) {
      m = ms[i];
      maxdet = det;
    }
  }
  return m;
}

class Mesh;

Int get_metric_dim(Int ncomps);
Int get_metrics_dim(LO nmetrics, Reals metrics);
Int get_metric_dim(Mesh* mesh);

Reals get_mident_metrics(
    Mesh* mesh, Int ent_dim, LOs entities, Reals v2m, bool has_degen = false);
Reals get_mident_metrics(
    Mesh* mesh, Int ent_dim, Reals v2m, bool has_degen = false);
Reals interpolate_between_metrics(LO nmetrics, Reals a, Reals b, Real t);
Reals linearize_metrics(LO nmetrics, Reals metrics);
Reals delinearize_metrics(LO nmetrics, Reals linear_metrics);

Reals project_metrics(Mesh* mesh, Reals e2m);

Reals clamp_metrics(LO nmetrics, Reals metrics, Real h_min, Real h_max);
Reals get_pure_implied_isos(Mesh* mesh);
Reals get_implied_isos(Mesh* mesh);
Reals get_element_implied_length_metrics(Mesh* mesh);
Reals get_pure_implied_metrics(Mesh* mesh);
Reals get_implied_metrics(Mesh* mesh);
Reals limit_metric_gradation(Mesh* mesh, Reals values, Real max_rate,
    Real tol = 1e-2, bool verbose = true);
Reals get_complexity_per_elem(Mesh* mesh, Reals v2m);
Reals get_nelems_per_elem(Mesh* mesh, Reals v2m);
Real get_complexity(Mesh* mesh, Reals v2m);
Real get_expected_nelems_from_complexity(Real complexity, Int dim);
Real get_expected_nelems(Mesh* mesh, Reals v2m);
Real get_metric_scalar_for_nelems(
    Int elem_dim, Real expected_nelems, Real target_nelems);
Real get_metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems);
Reals smooth_metric_once(Mesh* mesh, Reals v2m, bool has_dege = false);
Reals get_curvature_metrics(Mesh* mesh, Real segment_angle);
Reals get_hessian_metrics(Int dim, Reals hessians, Real eps);
Reals get_gradient_metrics(Int dim, Reals gradients, Real eps);
Reals get_derivative_metrics(Mesh* mesh, std::string const& name, Real knob);
Reals get_variation_metrics(Mesh* mesh, std::string const& name, Real knob);
Reals intersect_metrics(LO nmetrics, Reals a, Reals b);
Reals metrics_from_isos(Int new_dim, Reals isos);
Reals apply_isotropy(LO nmetrics, Reals metrics, Omega_h_Isotropy isotropy);

Reals isos_from_lengths(Reals h);
Reals lengths_from_isos(Reals l);

}  // end namespace Omega_h

#endif
