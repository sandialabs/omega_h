#ifndef FIT_HPP
#define FIT_HPP

/* This code is used to fit a linear polynomial to
 * element-wise sample points in a cavity.
 *
 * We do this by solving a least-squares problem using
 * a QR decomposition.
 *
 * MaxFitPoints denotes the maximum number of adjacent
 * sample points that will be accepted into the fitting
 * process.
 * These limits were chosen to be approximately twice the
 * number of elements adjacent to a vertex.
 */

#include "Omega_h_adj.hpp"
#include "Omega_h_qr.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

template <Int dim>
struct MaxFitPoints {
  enum { value = SimplexAvgDegree<dim, 0, dim>::value * 3 };
};

/* Computes the QR decomposition for the Vandermonde
 * matrix involved in the least-squares fit.
 * This depends only on the centroids of the adjacent
 * elements, nothing else.
 */

template <Int dim>
OMEGA_H_DEVICE QRFactorization<MaxFitPoints<dim>::value, dim + 1>
get_cavity_qr_factorization(LO k, LOs const& k2ke, LOs const& ke2e,
    LOs const& ev2v, Reals const& coords) {
  constexpr auto max_fit_pts = MaxFitPoints<dim>::value;
  Matrix<max_fit_pts, dim + 1> vandermonde;
  auto begin = k2ke[k];
  auto end = k2ke[k + 1];
  auto nfit_pts = end - begin;
  OMEGA_H_CHECK(nfit_pts >= dim + 1);
  OMEGA_H_CHECK(nfit_pts <= max_fit_pts);
  for (auto i = 0; i < nfit_pts; ++i) {
    auto ke = begin + i;
    auto e = ke2e[ke];
    auto eev2v = gather_verts<dim + 1>(ev2v, e);
    auto eev2x = gather_vectors<dim + 1, dim>(coords, eev2v);
    auto centroid = average(eev2x);
    vandermonde[0][i] = 1.0;
    for (Int j = 0; j < dim; ++j) {
      vandermonde[1 + j][i] = centroid[j];
    }
  }
  auto qr = factorize_qr_householder(nfit_pts, dim + 1, vandermonde);
// workaround CUDA compiler bug that makes these NaN
// if and only if we do not check whether they are NaN...
#ifdef __CUDA_ARCH__
  for (int k = 0; k < dim + 1; ++k) {
    for (int i = k; i < nfit_pts; ++i) {
      assert(!isnan(qr.v[k][i]));
    }
  }
#endif
  return qr;
}

/* Computes the linear polynomial coefficients
 * to approximat the distribution of one scalar
 * value over the cavity.
 * Typically this scalar is one of several stored
 * contiguously, hence the "ncomps" and "comp" arguments.
 * "e_data" contains the scalar data for all elements.
 * Notice that this function re-uses a pre-computed QR decomposition.
 */

template <Int dim>
OMEGA_H_DEVICE Vector<dim + 1> fit_cavity_polynomial(
    QRFactorization<MaxFitPoints<dim>::value, dim + 1> qr, LO k,
    LOs const& k2ke, LOs const& ke2e, Reals const& e_data, Int comp,
    Int ncomps) {
  constexpr auto max_fit_pts = MaxFitPoints<dim>::value;
  auto begin = k2ke[k];
  auto end = k2ke[k + 1];
  auto nfit_pts = end - begin;
  OMEGA_H_CHECK(nfit_pts >= dim + 1);
  OMEGA_H_CHECK(nfit_pts <= max_fit_pts);
  Vector<max_fit_pts> b;
  for (auto i = 0; i < nfit_pts; ++i) {
    auto ke = i + begin;
    auto e = ke2e[ke];
    b[i] = e_data[e * ncomps + comp];
  }
  for (auto i = nfit_pts; i < max_fit_pts; ++i) {
    b[i] = 0.0;
  }
  auto qtb = implicit_q_trans_b(nfit_pts, dim + 1, qr.v, b);
  auto coeffs = solve_upper_triangular(dim + 1, qr.r, qtb);
  return coeffs;
}

template <Int dim>
OMEGA_H_DEVICE Real eval_polynomial(Vector<dim + 1> coeffs, Vector<dim> x) {
  auto val = coeffs[0];
  for (Int j = 0; j < dim; ++j) val += coeffs[1 + j] * x[j];
  return val;
}

}  // end namespace Omega_h

#endif
