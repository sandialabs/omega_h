#ifndef METRIC_HPP
#define METRIC_HPP

#include "eigen.hpp"
#include "space.hpp"

namespace Omega_h {

template <Int dim>
INLINE Real metric_product(Matrix<dim, dim> m, Vector<dim> v) {
  return v * (m * v);
}

template <Int dim>
INLINE Real metric_length(Matrix<dim, dim> m, Vector<dim> v) {
  return sqrt(metric_product(m, v));
}

template <Int dim>
INLINE Real metric_desired_length(Matrix<dim, dim> m, Vector<dim> dir) {
  return 1.0 / metric_length(m, dir);
}

template <Int dim>
INLINE Vector<dim> metric_lengths(Vector<dim> l) {
  Vector<dim> h;
  for (Int i = 0; i < dim; ++i) h[i] = 1.0 / sqrt(l[i]);
  return h;
}

template <Int dim>
INLINE Decomposition<dim> decompose_metric(Matrix<dim, dim> m) {
  auto ed = decompose_eigen(m);
  auto h = metric_lengths(ed.l);
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
INLINE Matrix<dim, dim> intersect_metrics(
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
    Real u = metric_product(m1, p[i]);
    Real v = metric_product(m2, p[i]);
    w[i] = max2(u, v);
  }
  auto ip = invert(p);
  auto m = transpose(ip) * diagonal(w) * ip;
  return m;
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
That leaves (3) as being the best choice for these three reasons:
 - It has decent output in anisotropic cases
 - It can be used in triangles and tets
 - It does not require an eigendecomposition

Looking a (1), (2) and (3) suggests that their only
difference is an operation we will call "linearization",
in this case converting the metric tensor into a quantity
that can be safely linearly interpolated.
(1) M^{-1/2}
(2) M^{-1}
(3) M
*/

template <Int dim>
INLINE Matrix<dim, dim> linearize_metric(Matrix<dim, dim> m) {
  return log(m);
}

template <Int dim>
INLINE Matrix<dim, dim> delinearize_metric(Matrix<dim, dim> log_m) {
  return exp(log_m);
}

INLINE Real linearize_metric(Real h) { return ::log(h); }

INLINE Real delinearize_metric(Real log_h) { return ::exp(log_h); }

template <typename T>
INLINE T interpolate_metric(T a, T b, Real t) {
  return delinearize_metric(
      (linearize_metric(a) * (1.0 - t)) + (linearize_metric(b) * t));
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
INLINE Matrix<dim, dim> maxdet_metric(Few<Matrix<dim, dim>, n> ms) {
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

Reals get_midedge_metrics(Mesh* mesh, LOs entities, Reals v2m);
Reals get_maxdet_metrics(Mesh* mesh, LOs entities, Reals v2m);
Reals interpolate_between_metrics(Int dim, Reals a, Reals b, Real t);
Reals linearize_metrics(Int dim, Reals metrics);
Reals delinearize_metrics(Int dim, Reals linear_metrics);

Reals metric_from_hessians(
    Int dim, Reals hessians, Real eps, Real hmin, Real hmax);
Reals metric_for_nelems_from_hessians(Mesh* mesh, Real target_nelems,
    Real tolerance, Reals hessians, Real hmin, Real hmax);

/* used to achieve templated versions of code that either
 * accepts a metric tensor or nothing (nothing being the case
 * of isotropic quality, where the actual isotropic value doesn't
 * matter because our shape measure is scale-invariant)
 */
struct DummyIsoMetric {};

}  // end namespace Omega_h

#endif
