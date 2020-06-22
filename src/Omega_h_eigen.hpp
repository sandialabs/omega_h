#ifndef OMEGA_H_EIGEN_HPP
#define OMEGA_H_EIGEN_HPP

#include <Omega_h_matrix.hpp>

namespace Omega_h {

template <Int degree>
struct Roots {
  Int n;  // number of unique roots
  Few<Real, degree> values;
  Few<Int, degree> mults;  // multiplicities
};

// solve cubic equation x^3 + a_2 * x^2 + a_1 * x + a_0 = 0
// this code assumes that the solution does not have complex roots !
// the return value is the number of distinct roots,
// the two output arrays contain root values and multiplicities.
// roots within an *absolute* distance of (eps) are considered
// the same.
OMEGA_H_INLINE Roots<3> find_polynomial_roots(
    Few<Real, 3> coeffs, Real eps = 1e-6) {
  auto a_0 = coeffs[0];
  auto a_1 = coeffs[1];
  auto a_2 = coeffs[2];
  Few<Real, 3> roots;
  roots[0] = roots[1] = roots[2] = 0.0;
  Few<Int, 3> mults;
  mults[0] = mults[1] = mults[2] = 0;
  // http://mathworld.wolfram.com/CubicFormula.html
  Real p = (3. * a_1 - square(a_2)) / 3.;
  Real q = (9. * a_1 * a_2 - 27. * a_0 - 2. * cube(a_2)) / 27.;
  Real Q = p / 3.;
  Real R = q / 2.;
  Real D = cube(Q) + square(R);
  Real shift = -a_2 / 3.;
  if (D >= 0.0) {
    Real S = std::cbrt(R + std::sqrt(D));
    Real T = std::cbrt(R - std::sqrt(D));
    Real B = S + T;
    Real z_1 = shift + B;
    Real z_23_real = shift - (1. / 2.) * B;
    // warning: we simply assume the imaginary component is small !
    // Real z_23_imag = (1. / 2.) * sqrt(3.) * A;
    roots[0] = z_1;
    roots[1] = roots[2] = z_23_real;
  } else {
    // D < 0 implies Q < 0, since R^2 must be positive
    auto cos_theta = R / std::sqrt(-cube(Q));
    Real theta = std::acos(clamp(cos_theta, -1.0, 1.0));
    Real radius = 2. * std::sqrt(-Q);
    Real z_1 = radius * std::cos((theta) / 3.) + shift;
    Real z_2 = radius * std::cos((theta + 2. * PI) / 3.) + shift;
    Real z_3 = radius * std::cos((theta - 2. * PI) / 3.) + shift;
    roots[0] = z_1;
    roots[1] = z_2;
    roots[2] = z_3;
  }
  mults[0] = mults[1] = mults[2] = 1;
  // post-processing, decide whether pairs of roots are
  // close enough to be called repeated roots
  // first step, if two roots are close, then
  // move them to the second and third slots
  if (std::abs(roots[0] - roots[1]) < eps) {
    swap2(roots[0], roots[2]);
  } else if (std::abs(roots[0] - roots[2]) < eps) {
    swap2(roots[0], roots[1]);
  } else if (std::abs(roots[1] - roots[2]) < eps) {
    // no need to swap, they're already there
  } else {
    // no pairs were close, all three roots are distinct
    return {3, roots, mults};
  }
  // if we're here, two close roots are in [1] and [2]
  roots[1] = average(roots[1], roots[2]);
  mults[1] = 2;
  // lets see if they are all the same
  if (std::abs(roots[0] - roots[1]) < eps) {
    // roots[1] is already an average, weight it properly
    roots[0] = (1. / 3.) * roots[0] + (2. / 3.) * roots[1];
    mults[0] = 3;
    return {1, roots, mults};
  }
  return {2, roots, mults};
}

// solve quadratic equation x^2 + a * x + b = 0
OMEGA_H_INLINE Roots<2> find_polynomial_roots(
    Few<Real, 2> coeffs, Real eps = 1e-6) {
  auto a = coeffs[1];
  auto b = coeffs[0];
  Few<Real, 2> roots;
  roots[0] = roots[1] = 0.0;
  Few<Int, 2> mults;
  mults[0] = mults[1] = 0;
  Real disc = square(a) - 4. * b;
  if (std::abs(disc) < eps) {
    mults[0] = 2;
    roots[0] = -a / 2.;
    roots[1] = roots[0];
    mults[1] = mults[0];
    return {1, roots, mults};
  }
  if (disc > 0.0) {
    mults[0] = 1;
    mults[1] = 1;
    roots[0] = (-a + std::sqrt(disc)) / 2.;
    roots[1] = (-a - std::sqrt(disc)) / 2.;
    return {2, roots, mults};
  }
  return {0, roots, mults};
}

// solve linear equation x + a = 0
OMEGA_H_INLINE Roots<1> find_polynomial_roots(
    Few<Real, 1> coeffs, Real eps = 1e-6) {
  (void)eps;
  auto a = coeffs[0];
  Few<Real, 1> roots;
  roots[0] = -a;
  Few<Int, 1> mults;
  mults[0] = 1;
  return {1, roots, mults};
}

/* http://mathworld.wolfram.com/CharacteristicPolynomial.html */
OMEGA_H_INLINE Few<Real, 3> characteristic_polynomial(Tensor<3> A) {
  auto tA = trace(A);
  Few<Real, 3> coeffs;
  coeffs[2] = -tA;
  coeffs[1] = (1. / 2.) * ((tA * tA) - trace(A * A));
  coeffs[0] = -determinant(A);
  return coeffs;
}

/* http://mathworld.wolfram.com/CharacteristicPolynomial.html */
OMEGA_H_INLINE Few<Real, 2> characteristic_polynomial(Tensor<2> A) {
  Few<Real, 2> coeffs;
  coeffs[1] = -trace(A);
  coeffs[0] = determinant(A);
  return coeffs;
}

OMEGA_H_INLINE Few<Real, 1> characteristic_polynomial(Tensor<1> A) {
  Few<Real, 1> coeffs;
  coeffs[0] = -determinant(A);
  return coeffs;
}

template <Int n>
OMEGA_H_INLINE Roots<n> get_eigenvalues(Tensor<n> A) {
  auto poly = characteristic_polynomial(A);
  // WARNING: I no longer remember the source of this magic number.
  // was probably tuned to avoid failures with the cubic solver
  return find_polynomial_roots(poly, 5e-5);
}

template <Int m>
OMEGA_H_INLINE Tensor<m> subtract_from_diag(Tensor<m> a, Real const mu) {
  for (Int i = 0; i < m; ++i) a[i][i] -= mu;
  return a;
}

/* the null space of the matrix (s = m - l*I)
   is the space spanned by the eigenvectors.
   the multiplicity of this root is the dimensionality
   of that space, i.e. the number of eigenvectors.
   since the null space is the orthogonal component of
   the row span, we essentially first find the row space */

/* in the case that the null space is 1D and space is 3D,
   then the row space is 2D (a plane in 3D)
   take the largest cross product of any pair of rows,
   that should give a vector in the null space */
OMEGA_H_INLINE Vector<3> single_eigenvector(Tensor<3> const m, Real const l) {
  auto s = transpose(subtract_from_diag(m, l));
  auto v = cross(s[0], s[1]);
  auto v_norm = norm(v);
  auto c = cross(s[1], s[2]);
  auto c_norm = norm(c);
  if (c_norm > v_norm) {
    v = c;
    v_norm = c_norm;
  }
  c = cross(s[0], s[2]);
  c_norm = norm(c);
  if (c_norm > v_norm) {
    v = c;
    v_norm = c_norm;
  }
  OMEGA_H_CHECK(v_norm > EPSILON);
  v = v / v_norm;
  return v;
}

/* in the case that all rows of (a) are linearly dependent,
   this function will return the unit vector of the row
   that had the highest norm, which is a basis for the row space */
template <Int m>
OMEGA_H_INLINE Vector<m> get_1d_row_space(Tensor<m> const a) {
  auto ta = transpose(a);
  auto best_row = 0;
  auto best_norm = norm(ta[best_row]);
  for (Int i = 1; i < m; ++i) {
    auto row_norm = norm(ta[i]);
    if (row_norm > best_norm) {
      best_row = i;
      best_norm = row_norm;
    }
  }
  OMEGA_H_CHECK(best_norm > EPSILON);
  return ta[best_row] / best_norm;
}

/* in the case that the null space is 2D and space is 3D,
   find two vectors that are orthogonal to the 1D row space */
OMEGA_H_INLINE Few<Vector<3>, 2> double_eigenvector(
    Tensor<3> const m, Real const l) {
  auto s = subtract_from_diag(m, l);
  auto n = get_1d_row_space(s);
  auto b = form_ortho_basis(n);
  Few<Vector<3>, 2> o;
  o[0] = b[1];
  o[1] = b[2];
  return o;
}

template <Int dim>
struct DiagDecomp {
  Tensor<dim> q;
  Vector<dim> l;
};

OMEGA_H_INLINE DiagDecomp<3> decompose_eigen_dim(Tensor<3> const m) {
  auto roots_obj = get_eigenvalues(m);
  auto nroots = roots_obj.n;
  auto roots = roots_obj.values;
  auto mults = roots_obj.mults;
  /* there are only a few output cases, see solve_cubic() */
  Tensor<3> q;
  Vector<3> l;
  if (nroots == 3) {
    for (Int i = 0; i < 3; ++i) {
      q[i] = single_eigenvector(m, roots[i]);
      l[i] = roots[i];
    }
  } else if (nroots == 2 && mults[1] == 2) {
    q[0] = single_eigenvector(m, roots[0]);
    l[0] = roots[0];
    auto dev = double_eigenvector(m, roots[1]);
    q[1] = dev[0];
    q[2] = dev[1];
    l[1] = l[2] = roots[1];
  } else {
    OMEGA_H_CHECK(nroots == 1 && mults[0] == 3);
    l[0] = l[1] = l[2] = roots[0];
    q = identity_matrix<3, 3>();
  }
  return {q, l};
}

/* in the case that the null space is 1D and space is 2D,
   find the basis vector for the 1D row space and rotate it 90 deg */
OMEGA_H_INLINE Vector<2> single_eigenvector(Tensor<2> const m, Real const l) {
  return perp(get_1d_row_space(subtract_from_diag(m, l)));
}

OMEGA_H_INLINE DiagDecomp<2> decompose_eigen_dim(Tensor<2> const m) {
  auto roots_obj = get_eigenvalues(m);
  auto nroots = roots_obj.n;
  auto roots = roots_obj.values;
  /* there are only a few output cases, see solve_quadratic() */
  Tensor<2> q;
  Vector<2> l;
  if (nroots == 2) {
    for (Int i = 0; i < 2; ++i) {
      q[i] = single_eigenvector(m, roots[i]);
      l[i] = roots[i];
    }
  } else {
    OMEGA_H_CHECK(nroots == 1 && roots_obj.mults[0] == 2);
    l[0] = l[1] = roots[0];
    q = identity_matrix<2, 2>();
  }
  return {q, l};
}

OMEGA_H_INLINE DiagDecomp<1> decompose_eigen_dim(Tensor<1> const m) {
  auto roots_obj = get_eigenvalues(m);
  auto roots = roots_obj.values;
  OMEGA_H_CHECK(are_close(roots[0], m[0][0]));
  return {matrix_1x1(1.0), vector_1(roots[0])};
}

/* decompose an m x m matrix (where m <= 3) into
   eigenvalues and eigenvectors.

   note that (q) in the output is a matrix whose
   columns are the right eigenvectors.
   hence it actually corresponds to (Q^{-T})
   in the eigendecomposition
     M = Q \Lambda Q^{-1}
   where Q is the change of basis matrix
   the output should satisfy
     m ~= transpose(q * diagonal(l) * invert(q)) */
template <Int dim>
OMEGA_H_INLINE DiagDecomp<dim> decompose_eigen(Tensor<dim> m) {
  /* the cubic solver is especially sensitive to dynamic
     range. what we can do is to normalize the input matrix
     and then re-apply that norm to the resulting roots */
  Real nm = max_norm(m);
  if (nm <= EPSILON) {
    /* this is the zero matrix... */
    return {identity_matrix<dim, dim>(), zero_vector<dim>()};
  }
  m = m / nm;
  auto decomp = decompose_eigen_dim(m);
  return {decomp.q, decomp.l * nm};
}

/* Q, again, being the matrix whose columns
   are the right eigenvectors, but not necessarily unitary */
template <Int dim>
OMEGA_H_INLINE Tensor<dim> compose_eigen(
    Tensor<dim> const q, Vector<dim> const l) {
  return q * diagonal(l) * invert(q);
}

/* like the above, but knowing Q is unitary,
   so the transpose is the inverse */
template <Int dim>
OMEGA_H_INLINE Tensor<dim> compose_ortho(Tensor<dim> const q, Vector<dim> l) {
  return q * diagonal(l) * transpose(q);
}

// R^N off-diagonal norm. Useful for SVD and other algorithms
// that rely on Jacobi-type procedures.
// \param a
// \return \f$ \sqrt(\sum_i \sum_{j, j\neq i} a_{ij}^2) \f$
template <Int dim>
OMEGA_H_INLINE Real norm_off_diag(Tensor<dim> const a) {
  Real s = 0.0;
  for (Int j = 0; j < dim; ++j) {
    for (Int i = 0; i < dim; ++i) {
      if (i != j) {
        s += square(a(i, j));
      }
    }
  }
  return std::sqrt(s);
}

// R^N arg max off-diagonal. Useful for SVD and other algorithms
// that rely on Jacobi-type procedures.
// \param a
// \return \f$ (p,q) = arg max_{i \neq j} |a_{ij}| \f$
template <Int dim>
OMEGA_H_INLINE Few<Int, 2> arg_max_off_diag(Tensor<dim> const a) {
  Int p = 0;
  Int q = 0;
  auto s = -1.0;
  for (Int j = 0; j < dim; ++j) {
    for (Int i = 0; i < dim; ++i) {
      auto s2 = std::abs(a(i, j));
      if (i != j && s2 > s) {
        p = i;
        q = j;
        s = s2;
      }
    }
  }
  Few<Int, 2> out;
  out[0] = min2(p, q);
  out[1] = max2(p, q);
  return out;
}

// Symmetric Schur algorithm for R^2.
// \param \f$ A = [f, g; g, h] \in S(2) \f$
// \return \f$ c, s \rightarrow [c, -s; s, c]\f diagonalizes A$
OMEGA_H_INLINE Vector<2> schur_sym(Real f, Real g, Real h) {
  Real c = 1.0;
  Real s = 0.0;
  if (g != 0.0) {
    Real t = (h - f) / (2.0 * g);
    if (t >= 0.0) {
      t = 1.0 / (std::sqrt(1.0 + square(t)) + t);
    } else {
      t = -1.0 / (std::sqrt(1.0 + square(t)) - t);
    }
    c = 1.0 / std::sqrt(1.0 + square(t));
    s = t * c;
  }
  return vector_2(c, s);
}

/* Apply Givens-Jacobi rotation on the left */
template <Int dim>
OMEGA_H_INLINE Tensor<dim> givens_left(
    Real const c, Real const s, Int const i, Int const k, Tensor<dim> a) {
  for (Int j = 0; j < dim; ++j) {
    auto t1 = a(i, j);
    auto t2 = a(k, j);
    a(i, j) = c * t1 - s * t2;
    a(k, j) = s * t1 + c * t2;
  }
  return a;
}

/* Apply Givens-Jacobi rotation on the right */
template <Int dim>
OMEGA_H_INLINE Tensor<dim> givens_right(
    Real const c, Real const s, Int const i, Int const k, Tensor<dim> a) {
  for (Int j = 0; j < dim; ++j) {
    auto t1 = a(j, i);
    auto t2 = a(j, k);
    a(j, i) = c * t1 - s * t2;
    a(j, k) = s * t1 + c * t2;
  }
  return a;
}

/* the following is the Classic Jacobi algorithm, copied
   from the MiniTensor package, which in turn is based on
   algorithm 8.4.2 in Matrix Computations, Golub & Van Loan 1996 */
template <Int dim>
OMEGA_H_INLINE DiagDecomp<dim> decompose_eigen_jacobi(
    Tensor<dim> a, Real const eps = DBL_EPSILON, Int max_iter = -1) {
  // Estimate based on random generation and linear regression.
  // Golub & Van Loan p 429 expect ~ dimension * log(dimension)
  if (max_iter == -1) max_iter = (5 * dim * dim) / 2;
  auto v = identity_matrix<dim, dim>();
  auto tol = eps * norm(a);
  Int iter = 0;
  while (norm_off_diag(a) > tol && iter < max_iter) {
    auto pq = arg_max_off_diag(a);
    auto p = pq[0];
    auto q = pq[1];
    auto f = a(p, p);
    auto g = a(p, q);
    auto h = a(q, q);
    auto cs = schur_sym(f, g, h);
    auto c = cs[0];
    auto s = cs[1];
    a = givens_left(c, s, p, q, a);
    a = givens_right(c, s, p, q, a);
    v = givens_right(c, s, p, q, v);
    ++iter;
  }
  return {v, diagonal(a)};
}

template <Int dim>
OMEGA_H_INLINE DiagDecomp<dim> sort_by_magnitude(DiagDecomp<dim> const dd) {
  Few<Int, dim> perm;
  for (Int i = 0; i < dim; ++i) perm[i] = i;
  for (Int i = 0; i < dim; ++i) {
    auto max_j = i;
    auto max_m = std::abs(dd.l[perm[max_j]]);
    for (Int j = i + 1; j < dim; ++j) {
      auto m = std::abs(dd.l[perm[j]]);
      if (m > max_m) max_j = j;
    }
    swap2(perm[i], perm[max_j]);
  }
  DiagDecomp<dim> out;
  for (Int i = 0; i < dim; ++i) out.l[i] = dd.l[perm[i]];
  for (Int i = 0; i < dim; ++i) out.q[i] = dd.q[perm[i]];
  return out;
}

// logarithm of a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> log_spd_old(
    Tensor<dim> const m) OMEGA_H_NOEXCEPT {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::log(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// exponential resulting in a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> exp_spd_old(
    Tensor<dim> const m) OMEGA_H_NOEXCEPT {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::exp(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// logarithm of a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> log_spd(Tensor<dim> const m) OMEGA_H_NOEXCEPT {
  auto decomp = decompose_eigen_jacobi(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::log(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// exponential resulting in a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> exp_spd(Tensor<dim> const m) OMEGA_H_NOEXCEPT {
  auto decomp = decompose_eigen_jacobi(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::exp(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// exponential resulting in a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE_BIG Tensor<dim> sqrt_spd(Tensor<dim> const m) OMEGA_H_NOEXCEPT {
  auto decomp = decompose_eigen_jacobi(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::sqrt(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

Reals get_max_eigenvalues(Int dim, Reals symms);

}  // end namespace Omega_h

#endif
