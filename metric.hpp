#ifndef METRIC_HPP
#define METRIC_HPP

namespace osh {

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
INLINE Vector<dim> metric_eigenvalues(Vector<dim> h) {
  Vector<dim> l;
  for (Int i = 0; i < dim; ++i) l[i] = 1.0 / square(h[i]);
  return l;
}

template <Int dim>
INLINE Vector<dim> metric_lengths(Vector<dim> l) {
  Vector<dim> h;
  for (Int i = 0; i < dim; ++i) h[i] = 1.0 / sqrt(l[i]);
  return h;
}

template <Int dim>
INLINE Matrix<dim, dim> compose_metric(Matrix<dim, dim> r, Vector<dim> h) {
  auto l = metric_eigenvalues(h);
  return r * diagonal(l) * transpose(r);
}

template <Int dim>
INLINE void decompose_metric(Matrix<dim, dim> m, Matrix<dim, dim>& r,
                             Vector<dim>& h) {
  Vector<dim> l;
  decompose_eigen(m, r, l);
  h = metric_lengths(l);
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
*/

template <Int dim>
INLINE Matrix<dim, dim> common_metric_basis(Matrix<dim, dim> a,
                                            Matrix<dim, dim> b) {
  auto c = invert(a) * b;
  Matrix<dim, dim> p;
  Vector<dim> l;
  decompose_eigen(c, p, l);
  return p;
}

template <Int dim>
INLINE Matrix<dim, dim> intersect_metrics(Matrix<dim, dim> a,
                                          Matrix<dim, dim> b) {
  auto p = common_metric_basis(a, b);
  Vector<dim> w;
  for (Int i = 0; i < dim; ++i) {
    Real u = metric_product(a, p[i]);
    Real v = metric_product(b, p[i]);
    w[i] = max2(u, v);
  }
  auto mi = compose_eigen(p, w);
  return mi;
}

/* Alauzet details four different ways to interpolate
   the metric tensor:

1) M(t) = ((1-t)M_1^{-1/2} + t M_2^{-1/2)^2

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
*/

template <Int dim>
INLINE Matrix<dim, dim> interpolate_metrics(Matrix<dim, dim> a,
                                            Matrix<dim, dim> b, Real t) {
  return invert((invert(a) * (1.0 - t)) + (invert(b) * t));
}

/* currently we intend only to interpolate metrics to the
   barycenter of an entity, hence this simplified code */
template <Int dim, Int n>
INLINE Matrix<dim, dim> average_metrics(Few<Matrix<dim, dim>, n> m) {
  auto am = zero_matrix<dim, dim>();
  for (Int i = 0; i < n; ++i) {
    am = am + invert(m[i]);
  }
  am = am / double(n);
  am = invert(am);
  return am;
}

Reals average_metric(Mesh* mesh, Int ent_dim, LOs entities, Reals v2m);
Reals interpolate_metrics(Int dim, Reals a, Reals b, Real t);

} //end namespace osh

#endif
