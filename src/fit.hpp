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

#include "access.hpp"
#include "qr.hpp"

namespace osh {

template <Int dim>
struct MaxFitPoints;

template <>
struct MaxFitPoints<2> {
  enum { value = 12 };
};

template <>
struct MaxFitPoints<3> {
  enum { value = 48 };
};

template <Int dim>
struct CavityQrDecomposition {
  enum { max_fit_pts = MaxFitPoints<dim>::value };
  Few<Vector<max_fit_pts>, dim + 1> householder_vecs;
  Matrix<dim + 1, dim + 1> r;
};

/* Computes the QR decomposition for the Vandermonde
 * matrix involved in the least-squares fit.
 * This depends only on the centroids of the adjacent
 * elements, nothing else.
 */

template <Int dim>
DEVICE CavityQrDecomposition<dim> get_cavity_qr_decomposition(LO k,
    LOs const& k2ke, LOs const& ke2e, LOs const& ev2v, Reals const& coords) {
  constexpr auto max_fit_pts = MaxFitPoints<dim>::value;
  Matrix<max_fit_pts, dim + 1> vandermonde;
  CHECK(k2ke[k + 1] - k2ke[k] <= max_fit_pts);
  Int nfit_pts = 0;
  for (auto ke = k2ke[k]; ke < k2ke[k + 1]; ++ke) {
    auto e = ke2e[ke];
    auto eev2v = gather_verts<dim + 1>(ev2v, e);
    auto eev2x = gather_vectors<dim + 1, dim>(coords, eev2v);
    auto centroid = average(eev2x);
    vandermonde[0][nfit_pts] = 1.0;
    for (Int j = 0; j < dim; ++j) {
      vandermonde[1 + j][nfit_pts] = centroid[j];
    }
    ++nfit_pts;
  }
  for (auto i = nfit_pts; i < max_fit_pts; ++i) {
    for (Int j = 0; j < (dim + 1); ++j) {
      vandermonde[j][i] = 0.0;
    }
  }
  Few<Vector<max_fit_pts>, dim + 1> householder_vecs;
  auto r_full = vandermonde;
  auto rank = factorize_qr_householder(r_full, householder_vecs);
  CHECK(rank == dim + 1);
  auto r = reduced_r_from_full(r_full);
  return {householder_vecs, r};
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
DEVICE Vector<dim + 1> fit_cavity_polynomial(
    CavityQrDecomposition<dim> qr_decomp, LO k, LOs const& k2ke,
    LOs const& ke2e, Reals const& e_data, Int comp, Int ncomps) {
  constexpr auto max_fit_pts = MaxFitPoints<dim>::value;
  Vector<max_fit_pts> b;
  Int nfit_pts = 0;
  for (auto ke = k2ke[k]; ke < k2ke[k + 1]; ++ke) {
    auto e = ke2e[ke];
    b[nfit_pts] = e_data[e * ncomps + comp];
    ++nfit_pts;
  }
  for (auto i = nfit_pts; i < max_fit_pts; ++i) b[i] = 0;
  auto qtb_full = b;
  implicit_q_trans_b(qtb_full, qr_decomp.householder_vecs);
  Vector<dim + 1> qtb;
  for (Int i = 0; i < dim + 1; ++i) qtb[i] = qtb_full[i];
  auto coeffs = solve_upper_triangular(qr_decomp.r, qtb);
  return coeffs;
}

template <Int dim>
DEVICE Real eval_polynomial(Vector<dim + 1> coeffs, Vector<dim> x) {
  auto val = coeffs[0];
  for (Int j = 0; j < dim; ++j) val += coeffs[1 + j] * x[j];
  return val;
}

}  // end namespace osh

#endif
