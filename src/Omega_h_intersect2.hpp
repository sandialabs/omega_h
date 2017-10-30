#include <Omega_h_eigen.hpp>

OMEGA_H_INLINE Matrix<1, 1> intersect_degenerate_metrics(
    Matrix<1, 1> m1,
    DiagDecomp<1, 1>,
    Few<Int, 1>,
    Int,
    Matrix<1, 1>,
    DiagDecomp<1, 1>,
    Few<Int, 1>,
    Int) {
  // this should be impossible, but toss something in here so the code compiles
  return m1;
}

// Appendix A.1 in Barral's dissertation, case 5
OMEGA_H_INLINE Matrix<2, 2> intersect_degenerate_metrics(
    Matrix<2, 2> m1,
    DiagDecomp<2, 2> m1_dc,
    Few<Int, 2> m1_ew_is_degen,
    Int,
    Matrix<2, 2> m2,
    DiagDecomp<2, 2> m2_dc,
    Few<Int, 2> m2_ew_is_degen,
    Int) {
  Vector<2> u1, v1, v2;
  Real l1, l2;
  for (Int i = 0; i < 2; ++i) {
    if (m1_ew_is_degen[i]) v1 = m1_dc.q[i];
    else {
      u1 = m1_dc.q[i];
      l1 = m1_dc.l[i];
    }
    if (m2_ew_is_degen[i]) v2 = m2_dc.q[i];
    else l2 = m2_dc.l[i];
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
    l[0] = v1 * (M2 * v1);
    l[1] = v2 * (M1 * v2);
    return tranpose(p_inv) * diagonal(l) * p_inv;
  }
}

OMEGA_H_INLINE Matrix<3, 3> intersect_degenerate_metrics(
    Matrix<3, 3> m1,
    DiagDecomp<3, 3> m1_dc,
    Few<Int, 3> m1_ew_is_degen,
    Int nm1_degen_ews,
    Matrix<3, 3> m2,
    DiagDecomp<3, 3> m2_dc,
    Few<Int, 3> m2_ew_is_degen,
    Int nm2_degen_ews) {
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> intersect_metrics(
    Matrix<dim, dim> m1,
    Matrix<dim, dim> m2) {
  auto m1_dc = decompose_eigen(m1);
  auto m2_dc = decompose_eigen(m2);
  Few<Int, dim> m1_ew_is_degen;
  Few<Int, dim> m2_ew_is_degen;
  for (Int i = 0; i < dim; ++i) {
    // if EPSILON is 1e-10, this corresponds roughly to edges
    // of absolute length 1e5 or more being considered "infinite"
    m1_ew_is_degen = m1_dc.l[i] < OMEGA_H_EPSILON;
    m2_ew_is_degen = m2_dc.l[i] < OMEGA_H_EPSILON;
  }
  auto nm1_degen_ews = sum(m1_ew_is_degen);
  auto nm2_degen_ews = sum(m2_ew_is_degen);
  if (nm1_degen_ews > nm2_degen_ews) {
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
    auto m2_bar = tranpose(m1_inv_sqrt) * m2 * m1_inv_sqrt;
    auto m2_bar_dc = decompose_eigen(m2_bar);
    auto p = m2_bar_dc.q;
    Vector<dim> l_m_int_bar;
    for (Int i = 0; i < dim; ++i) {
      l_m_int_bar[i] = max2(m2_bar_dc.l[i], 1.0);
    }
    auto m_int_bar = compose_ortho(p, l_m_int_bar);
    auto m_int = tranpose(m1_sqrt) * m_int_bar * m1_sqrt;
    return m_int;
  }
  // okay, both the metrics are partially degenerate.
  // call the dimension-specific logic.
}
