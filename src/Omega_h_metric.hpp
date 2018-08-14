#ifndef OMEGA_H_METRIC_HPP
#define OMEGA_H_METRIC_HPP

#include <Omega_h_eigen.hpp>
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
  return std::sqrt(metric_product(m, v));
}

template <Int dim>
OMEGA_H_INLINE Real metric_desired_length(Matrix<dim, dim> m, Vector<dim> dir) {
  return 1.0 / metric_length(m, dir);
}

OMEGA_H_INLINE Real metric_length_from_eigenvalue(Real l) {
  return 1.0 / std::sqrt(l);
}

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

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> intersect_metrics(
    Matrix<dim, dim> m1, Matrix<dim, dim> m2);

/* Metric intersection that accounts for all degenerate cases,
   thanks to:
   Barral, Nicolas.
   "Time-accurate anisotropic mesh adaptation for three-dimensional moving mesh
   problems" Diss. Universite Pierre et Marie Curie-Paris VI, 2015.
 */
OMEGA_H_INLINE Matrix<1, 1> intersect_degenerate_metrics(Matrix<1, 1> m1,
    DiagDecomp<1>, Few<Int, 1>, Int, Matrix<1, 1>, DiagDecomp<1>, Few<Int, 1>,
    Int) {
  // this should be impossible, but toss something in here so the code compiles
  return m1;
}

// Appendix A.1 in Barral's dissertation, case 5
OMEGA_H_INLINE Matrix<2, 2> intersect_degenerate_metrics(Matrix<2, 2> m1,
    DiagDecomp<2> m1_dc, Few<Int, 2> m1_ew_is_degen, Int, Matrix<2, 2> m2,
    DiagDecomp<2> m2_dc, Few<Int, 2> m2_ew_is_degen, Int) {
  auto u1 = zero_vector<2>();
  auto v1 = zero_vector<2>();
  auto v2 = zero_vector<2>();
  Real l1 = -1.0;
  Real l2 = -1.0;
  for (Int i = 0; i < 2; ++i) {
    if (m1_ew_is_degen[i])
      v1 = m1_dc.q[i];
    else {
      u1 = m1_dc.q[i];
      l1 = m1_dc.l[i];
    }
    if (m2_ew_is_degen[i])
      v2 = m2_dc.q[i];
    else
      l2 = m2_dc.l[i];
  }
  if (v1 * v2 < OMEGA_H_EPSILON) {
    // case 5.a
    return outer_product(u1, u1) * max2(l1, l2);
  } else {
    // case 5.b
    Matrix<2, 2> p;
    p[0] = v1;
    p[1] = v2;
    auto p_inv = invert(p);
    Vector<2> l;
    l[0] = v1 * (m2 * v1);
    l[1] = v2 * (m1 * v2);
    return transpose(p_inv) * diagonal(l) * p_inv;
  }
}

// Barral's thesis, appendix A.2
OMEGA_H_INLINE Matrix<3, 3> intersect_degenerate_metrics(Matrix<3, 3> m1,
    DiagDecomp<3> m1_dc, Few<Int, 3> m1_ew_is_degen, Int nm1_degen_ews,
    Matrix<3, 3> m2, DiagDecomp<3> m2_dc, Few<Int, 3> m2_ew_is_degen,
    Int nm2_degen_ews) {
  if (nm1_degen_ews == 2 && nm2_degen_ews == 2) {
    // case 2
    auto u1 = zero_vector<3>();
    auto u2 = zero_vector<3>();
    Real l1 = -1.0, l2 = -1.0;
    for (Int i = 0; i < 3; ++i) {
      if (!m1_ew_is_degen[i]) {
        u1 = m1_dc.q[i];
        l1 = m1_dc.l[i];
      }
      if (!m2_ew_is_degen[i]) {
        u2 = m2_dc.q[i];
        l2 = m2_dc.l[i];
      }
    }
    auto u = cross(u1, u2);
    if (norm_squared(u) < EPSILON) {
      // case 2.a (u1 == u2)
      return max2(l1, l2) * outer_product(u1, u1);
    } else {
      u = normalize(u);
      // case 2.b (u1 != u2)
      auto e1 = cross(u1, u);
      auto e2 = cross(u2, u);
      Matrix<3, 3> P;
      P[0] = e1;
      P[1] = e2;
      P[2] = u;
      Vector<3> l;
      l[0] = e1 * (m2 * e1);
      l[1] = e2 * (m1 * e2);
      l[2] = 0.0;
      auto Pinv = invert(P);
      return transpose(Pinv) * diagonal(l) * Pinv;
    }
  }
  if (nm1_degen_ews == 1 && nm2_degen_ews == 2) {
    // case 3
    // note that in Barral's dissertation, it is m1 that has two degenerate
    // directions. however, here we keep the rule that m1 is the least
    // degenerate. so, in this case all 1 and 2 are swapped compared to the
    // dissertation.
    auto u1 = zero_vector<3>();
    auto v1 = zero_vector<3>();
    auto w1 = zero_vector<3>();
    auto u2 = zero_vector<3>();
    auto v2 = zero_vector<3>();
    auto w2 = zero_vector<3>();
    bool found_u1 = false;
    bool found_v2 = false;
    for (Int i = 0; i < 3; ++i) {
      if (m1_ew_is_degen[i])
        w1 = m1_dc.q[i];
      else {
        if (found_u1) {
          v1 = m1_dc.q[i];
        } else {
          u1 = m1_dc.q[i];
          found_u1 = true;
        }
      }
      if (!m2_ew_is_degen[i])
        u2 = m2_dc.q[i];
      else {
        if (found_v2) {
          w2 = m2_dc.q[i];
        } else {
          v2 = m2_dc.q[i];
          found_v2 = true;
        }
      }
    }
    if (u2 * w1 < EPSILON) {
      // case 3.a , u2 and w1 are orthogonal
      Matrix<3, 2> P;
      P[0] = u1;
      P[1] = v1;
      auto PT = transpose(P);
      auto m1_bar = PT * (m1 * P);
      auto m2_bar = PT * (m2 * P);
      auto mint_bar = intersect_metrics(m1_bar, m2_bar);  // reduced to 2D
      // u1 and v1 are supposed to be orthogonal, so the pseudo-inverse is the
      // transpose
      return P * (mint_bar * PT);
    } else {
      // case 3.b, u2 and w1 are not orthogonal
      Matrix<3, 3> P;
      P[0] = v2;
      P[1] = w2;
      P[2] = w1;
      auto Pinv = invert(P);
      Vector<3> l;
      l[0] = P[0] * (m1 * P[0]);
      l[1] = P[1] * (m1 * P[1]);
      l[2] = P[2] * (m2 * P[2]);
      return transpose(Pinv) * diagonal(l) * Pinv;
    }
  }
  if (nm1_degen_ews == 1 && nm2_degen_ews == 1) {
    // case 4
    auto u1 = zero_vector<3>();
    auto v1 = zero_vector<3>();
    auto w1 = zero_vector<3>();
    auto w2 = zero_vector<3>();
    bool found_u1 = false;
    for (Int i = 0; i < 3; ++i) {
      if (m1_ew_is_degen[i])
        w1 = m1_dc.q[i];
      else {
        if (found_u1) {
          v1 = m1_dc.q[i];
        } else {
          u1 = m1_dc.q[i];
          found_u1 = true;
        }
      }
      if (m2_ew_is_degen[i]) w2 = m2_dc.q[i];
    }
    auto w = cross(w1, w2);
    if (norm_squared(w) < EPSILON) {
      // case 4.a
      Matrix<3, 2> P;
      P[0] = u1;
      P[1] = v1;
      auto PT = transpose(P);
      auto m1_bar = PT * (m1 * P);
      auto m2_bar = PT * (m2 * P);
      auto mint_bar = intersect_metrics(m1_bar, m2_bar);
      return P * (mint_bar * PT);
    } else {
      // case 4.b
      w = normalize(w);
      Matrix<3, 3> P;
      P[0] = w1;
      P[1] = w2;
      P[2] = w;
      Vector<3> l;
      l[0] = P[0] * (m2 * P[0]);
      l[1] = P[1] * (m1 * P[1]);
      l[2] = max2(P[2] * (m1 * P[2]), P[2] * (m2 * P[2]));
      auto Pinv = invert(P);
      return transpose(Pinv) * diagonal(l) * Pinv;
    }
  }
  return m1;
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> intersect_metrics(
    Matrix<dim, dim> m1, Matrix<dim, dim> m2) {
  auto m1_dc = decompose_eigen(m1);
  auto m2_dc = decompose_eigen(m2);
  Few<Int, dim> m1_ew_is_degen;
  Few<Int, dim> m2_ew_is_degen;
  for (Int i = 0; i < dim; ++i) {
    // if EPSILON is 1e-10, this corresponds roughly to edges
    // of absolute length 1e5 or more being considered "infinite"
    m1_ew_is_degen[i] = m1_dc.l[i] < OMEGA_H_EPSILON;
    m2_ew_is_degen[i] = m2_dc.l[i] < OMEGA_H_EPSILON;
  }
  auto nm1_degen_ews = reduce(m1_ew_is_degen, plus<Int>());
  auto nm2_degen_ews = reduce(m2_ew_is_degen, plus<Int>());
  if (nm1_degen_ews > nm2_degen_ews) {
    swap2(m1, m2);
    swap2(m1_dc, m2_dc);
    swap2(m1_ew_is_degen, m2_ew_is_degen);
    swap2(nm1_degen_ews, nm2_degen_ews);
  }
  // At this point, M_1 is the least degenerate, or they are equally degenerate.
  if (nm1_degen_ews == dim) {
    // The least degenerate metric is null... they must both be... return null
    return zero_matrix<dim, dim>();
  }
  if (nm1_degen_ews == 0) {
    // as long as M_1 is not degenerate, the "normal" procedure can be applied
    Vector<dim> l_m1_sqrt;
    Vector<dim> l_m1_inv_sqrt;
    for (Int i = 0; i < dim; ++i) {
      l_m1_sqrt[i] = std::sqrt(m1_dc.l[i]);
      l_m1_inv_sqrt[i] = 1.0 / l_m1_sqrt[i];
    }
    auto m1_sqrt = compose_ortho(m1_dc.q, l_m1_sqrt);
    auto m1_inv_sqrt = compose_ortho(m1_dc.q, l_m1_inv_sqrt);
    auto m2_bar = transpose(m1_inv_sqrt) * m2 * m1_inv_sqrt;
    auto m2_bar_dc = decompose_eigen(m2_bar);
    auto p = m2_bar_dc.q;
    Vector<dim> l_m_int_bar;
    for (Int i = 0; i < dim; ++i) {
      l_m_int_bar[i] = max2(m2_bar_dc.l[i], 1.0);
    }
    auto m_int_bar = compose_ortho(p, l_m_int_bar);
    auto m_int = transpose(m1_sqrt) * m_int_bar * m1_sqrt;
    return m_int;
  }
  // okay, both the metrics are partially degenerate.
  // call the dimension-specific logic.
  return intersect_degenerate_metrics(m1, m1_dc, m1_ew_is_degen, nm1_degen_ews,
      m2, m2_dc, m2_ew_is_degen, nm2_degen_ews);
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
template <Int dim>
OMEGA_H_INLINE void average_metric_contrib(
    Matrix<dim, dim>& am, Int& n, Matrix<dim, dim> m, bool has_degen) {
  if (has_degen && max_norm(m) < OMEGA_H_EPSILON) return;
  am += linearize_metric(m);
  n++;
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> average_metric_finish(
    Matrix<dim, dim> am, Int n, bool has_degen) {
  if (has_degen && n == 0) return am;
  am /= n;
  return delinearize_metric(am);
}

template <Int dim, Int n>
OMEGA_H_INLINE Matrix<dim, dim> average_metric(
    Few<Matrix<dim, dim>, n> ms, bool has_degen) {
  auto am = zero_matrix<dim, dim>();
  Int ngood = 0;
  for (Int i = 0; i < n; ++i) {
    average_metric_contrib(am, ngood, ms[i], has_degen);
  }
  return average_metric_finish(am, ngood, has_degen);
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
void axes_from_metric_field(
    Mesh* mesh, std::string const& metric_name, std::string const& axis_prefix);
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
