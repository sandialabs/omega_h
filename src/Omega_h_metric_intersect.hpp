#ifndef OMEGA_H_METRIC_INTERSECT_HPP
#define OMEGA_H_METRIC_INTERSECT_HPP

#include <Omega_h_metric.hpp>

namespace Omega_h {

template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> intersect_metrics(
    Tensor<dim> const m1, Tensor<dim> const m2);

/* Metric intersection that accounts for all degenerate cases,
   thanks to:
   Barral, Nicolas.
   "Time-accurate anisotropic mesh adaptation for three-dimensional moving mesh
   problems" Diss. Universite Pierre et Marie Curie-Paris VI, 2015.
 */
OMEGA_H_INLINE Tensor<1> intersect_degenerate_metrics(Tensor<1> const m1,
    DiagDecomp<1> const, Few<Int, 1> const, Int const, Tensor<1> const,
    DiagDecomp<1> const, Few<Int, 1> const, Int const) {
  // this should be impossible, but toss something in here so the code compiles
  return m1;
}

// Appendix A.1 in Barral's dissertation, case 5
OMEGA_H_INLINE_BIG Tensor<2> intersect_degenerate_metrics(Tensor<2> const m1,
    DiagDecomp<2> const m1_dc, Few<Int, 2> const m1_ew_is_degen, Int const,
    Tensor<2> const m2, DiagDecomp<2> const m2_dc,
    Few<Int, 2> const m2_ew_is_degen, Int const) {
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
    Tensor<2> p;
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
OMEGA_H_INLINE_BIG Tensor<3> intersect_degenerate_metrics(Tensor<3> const m1,
    DiagDecomp<3> const m1_dc, Few<Int, 3> const m1_ew_is_degen,
    Int const nm1_degen_ews, Tensor<3> const m2, DiagDecomp<3> const m2_dc,
    Few<Int, 3> const m2_ew_is_degen, Int const nm2_degen_ews) {
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
      Tensor<3> P;
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
      Tensor<3> P;
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
      Tensor<3> P;
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
OMEGA_H_INLINE_BIG Tensor<dim> intersect_metrics(
    Tensor<dim> const m1_in, Tensor<dim> const m2_in) {
  Tensor<dim> m1 = m1_in;
  Tensor<dim> m2 = m2_in;
  static_assert(std::is_same<decltype(m1), Tensor<dim>>::value, "not const!");
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
    swap2<Tensor<dim>>(m1, m2);
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

}  // namespace Omega_h

#endif
