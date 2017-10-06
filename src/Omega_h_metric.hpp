#ifndef OMEGA_H_METRIC_HPP
#define OMEGA_H_METRIC_HPP

#include <Omega_h_lie.hpp>

namespace Omega_h {

template <Int dim>
OMEGA_H_INLINE Real metric_product(Matrix<dim, dim> m, Vector<dim> v) {
  return v * (m * v);
}

template <Int space_dim>
OMEGA_H_INLINE typename std::enable_if<(space_dim > 1), Real>::type
metric_product(Matrix<1, 1> m, Vector<space_dim> v) {
  return v * (m[0][0] * v);
}

template <Int metric_dim, Int space_dim>
OMEGA_H_INLINE Real metric_length(
    Matrix<metric_dim, metric_dim> m, Vector<space_dim> v) {
  return sqrt(metric_product(m, v));
}

template <Int dim>
OMEGA_H_INLINE Real metric_desired_length(Matrix<dim, dim> m, Vector<dim> dir) {
  return 1.0 / metric_length(m, dir);
}

OMEGA_H_INLINE Real metric_length_from_eigenvalue(Real l)
    __attribute__((const));
OMEGA_H_INLINE Real metric_length_from_eigenvalue(Real l) {
  return 1.0 / sqrt(l);
}

OMEGA_H_INLINE Real metric_eigenvalue_from_length(Real h)
    __attribute__((const));
OMEGA_H_INLINE Real metric_eigenvalue_from_length(Real h) {
  return 1.0 / square(h);
}

template <Int dim>
OMEGA_H_INLINE Vector<dim> metric_lengths_from_eigenvalues(Vector<dim> l) {
  Vector<dim> h;
  for (Int i = 0; i < dim; ++i) h[i] = metric_length_from_eigenvalue(l[i]);
  return h;
}

template <Int dim>
OMEGA_H_INLINE Vector<dim> metric_eigenvalues_from_lengths(Vector<dim> h) {
  Vector<dim> l;
  for (Int i = 0; i < dim; ++i) l[i] = metric_eigenvalue_from_length(h[i]);
  return l;
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> compose_metric(
    Matrix<dim, dim> r, Vector<dim> h) {
  auto l = metric_eigenvalues_from_lengths(h);
  return compose_ortho(r, l);
}

template <Int dim>
OMEGA_H_INLINE DiagDecomp<dim> decompose_metric(Matrix<dim, dim> m) {
  auto ed = decompose_eigen(m);
  auto h = metric_lengths_from_eigenvalues(ed.l);
  return {ed.q, h};
}

/* INRIA knows what they're doing with respect to the metric.
   See these papers:

   Frey, Pascal-Jean, and Frédéric Alauzet.
   "Anisotropic mesh adaptation for CFD computations."
   Computer methods in applied mechanics and engineering
   194.48 (2005): 5068-5082.

   F. Alauzet, P.J. Frey, Estimateur d'erreur geometrique
   et metriques anisotropes pour l'adaptation de maillage.
   Partie I: aspects theoriques,
   RR-4759, INRIA Rocquencourt, 2003.

https://www.rocq.inria.fr/gamma/Frederic.Alauzet/

   This metric interpolation code is from Section 2.2 of
   that report, with one slight correction to Remark 2.5,
   where we suspect that in the more general case, if
   all eigenvalues of N are above one or all are below one,
   then one metric encompasses the other.
*/

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> intersect_metrics(
    Matrix<dim, dim> m1, Matrix<dim, dim> m2) {
  auto n = invert(m1) * m2;
  auto n_decomp = decompose_eigen(n);
  bool all_above_one = true;
  bool all_below_one = true;
  for (Int i = 0; i < dim; ++i) {
    if (n_decomp.l[i] > 1) all_below_one = false;
    if (n_decomp.l[i] < 1) all_above_one = false;
  }
  if (all_below_one) return m1;
  if (all_above_one) return m2;
  auto p = n_decomp.q;
  Vector<dim> w;
  for (Int i = 0; i < dim; ++i) {
    auto u = metric_product(m1, p[i]);
    auto v = metric_product(m2, p[i]);
    w[i] = max2(u, v);
  }
  auto ip = invert(p);
  return transpose(ip) * diagonal(w) * ip;
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
OMEGA_H_INLINE Matrix<dim, dim> linearize_metric(Matrix<dim, dim> m) {
  return log_spd(m);
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> delinearize_metric(Matrix<dim, dim> log_m) {
  return exp_spd(log_m);
}

template <Int n, typename T>
OMEGA_H_INLINE Few<T, n> linearize_metrics(Few<T, n> ms) {
  Few<T, n> log_ms;
  for (Int i = 0; i < n; ++i) log_ms[i] = linearize_metric(ms[i]);
  return log_ms;
}

/* the "proper" way to interpolate the metric tensor to
 * the barycenter of a simplex; does several eigendecompositions
 */
template <Int n, typename T>
OMEGA_H_INLINE T average_metric(Few<T, n> ms) {
  return delinearize_metric(average(linearize_metrics(ms)));
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> clamp_metric(
    Matrix<dim, dim> m, Real h_min, Real h_max) {
  auto ed = decompose_eigen(m);
  auto l_max = metric_eigenvalue_from_length(h_min);
  auto l_min = metric_eigenvalue_from_length(h_max);
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
OMEGA_H_INLINE Matrix<dim, dim> maxdet_metric(Few<Matrix<dim, dim>, n> ms) {
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

Reals get_mident_metrics(Mesh* mesh, Int ent_dim, LOs entities, Reals v2m);
Reals get_mident_metrics(Mesh* mesh, Int ent_dim, Reals v2m);
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
void axes_from_metric_field(
    Mesh* mesh, std::string const& metric_name, std::string const& axis_prefix);
Reals limit_metric_gradation(Mesh* mesh, Reals values, Real max_rate,
    Real tol = 1e-3, bool verbose = true);
Reals get_expected_nelems_per_elem(Mesh* mesh, Reals v2m);
Real get_expected_nelems(Mesh* mesh, Reals v2m);
Real get_metric_scalar_for_nelems(
    Int elem_dim, Real expected_nelems, Real target_nelems);
Real get_metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems);
Reals smooth_metric_once(Mesh* mesh, Reals v2m);
Reals get_curvature_metrics(Mesh* mesh, Real segment_angle);
Reals get_hessian_metrics(Int dim, Reals hessians, Real eps);
Reals get_gradient_metrics(Int dim, Reals gradients, Real eps);
Reals get_derivative_metrics(Mesh* mesh, std::string const& name, Real knob);
Reals get_variation_metrics(Mesh* mesh, std::string const& name, Real knob);
Reals get_proximity_isos(Mesh* mesh, Real factor);
Reals intersect_metrics(LO nmetrics, Reals a, Reals b);
Reals metrics_from_isos(Int new_dim, Reals isos);
Reals apply_isotropy(LO nmetrics, Reals metrics, Omega_h_Isotropy isotropy);
Reals get_aniso_zz_metric(
    Mesh* mesh, Reals elem_gradients, Real error_bound, Real max_size);

Reals isos_from_lengths(Reals h);
Reals lengths_from_isos(Reals l);

}  // end namespace Omega_h

#endif
