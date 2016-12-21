#ifndef EIGEN_HPP
#define EIGEN_HPP

#include "polynomial.hpp"
#include "space.hpp"

namespace Omega_h {

/* http://mathworld.wolfram.com/CharacteristicPolynomial.html */
INLINE Few<Real, 3> characteristic_polynomial(Matrix<3, 3> A) {
  auto tA = trace(A);
  Few<Real, 3> coeffs;
  coeffs[2] = -tA;
  coeffs[1] = (1. / 2.) * ((tA * tA) - trace(A * A));
  coeffs[0] = -determinant(A);
  return coeffs;
}

/* http://mathworld.wolfram.com/CharacteristicPolynomial.html */
INLINE Few<Real, 2> characteristic_polynomial(Matrix<2, 2> A) {
  Few<Real, 2> coeffs;
  coeffs[1] = -trace(A);
  coeffs[0] = determinant(A);
  return coeffs;
}

template <Int n>
INLINE Roots<n> get_eigenvalues(Matrix<n, n> A) {
  auto poly = characteristic_polynomial(A);
  return find_polynomial_roots(poly);
}

/* the null space of the matrix (s = m - l*I)
   is the space spanned by the eigenvectors.
   the multiplicity of this root is the dimensionality
   of that space, i.e. the number of eigenvectors. */

/* in the case that the null space is 1D and space is 3D,
   take the largest cross product of any pair of columns */
INLINE Vector<3> single_eigenvector(Matrix<3, 3> m, Real l) {
  subtract_from_diag(m, l);
  auto v = cross(m[0], m[1]);
  Real v_norm = norm(v);
  Vector<3> c = cross(m[1], m[2]);
  Real c_norm = norm(c);
  if (c_norm > v_norm) {
    v = c;
    v_norm = c_norm;
  }
  c = cross(m[0], m[2]);
  c_norm = norm(c);
  if (c_norm > v_norm) {
    v = c;
    v_norm = c_norm;
  }
  CHECK(v_norm > EPSILON);
  v = v / v_norm;
  return v;
}

template <Int m>
INLINE Vector<m> get_1d_column_space(Matrix<m, m> a) {
  Vector<m> v = zero_vector<m>();
  Real v_norm = 0;
  for (Int j = 0; j < m; ++j) {
    Real c_norm = norm(a[j]);
    if (c_norm > v_norm) {
      v = a[j];
      v_norm = c_norm;
    }
  }
  CHECK(v_norm > EPSILON);
  return v / v_norm;
}

/* in the case that the null space is 2D, find the
   largest-norm column and get a couple vectors
   orthogonal to that */
INLINE Few<Vector<3>, 2> double_eigenvector(Matrix<3, 3> m, Real l) {
  subtract_from_diag(m, l);
  Vector<3> n = get_1d_column_space(m);
  Matrix<3, 3> b = form_ortho_basis(n);
  Few<Vector<3>, 2> out;
  out[0] = b[1];
  out[1] = b[2];
  return out;
}

template <Int dim>
struct DiagDecomp {
  Matrix<dim, dim> q;
  Vector<dim> l;
};

INLINE DiagDecomp<3> decompose_eigen_dim(Matrix<3, 3> m) {
  auto roots_obj = get_eigenvalues(m);
  auto nroots = roots_obj.n;
  auto roots = roots_obj.values;
  auto mults = roots_obj.mults;
  /* there are only a few output cases, see solve_cubic() */
  Matrix<3, 3> q;
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
    CHECK(nroots == 1 && mults[0] == 3);
    l[0] = l[1] = l[2] = roots[0];
    q = identity_matrix<3, 3>();
  }
  return {q, l};
}

/* in the case that the null space is 1D and space is 2D,
   get the largest column and rotate it 90 deg */
INLINE Vector<2> single_eigenvector(Matrix<2, 2> m, Real l) {
  Matrix<2, 2> s = (m - (l * identity_matrix<2, 2>()));
  return perp(get_1d_column_space(s));
}

INLINE DiagDecomp<2> decompose_eigen_dim(Matrix<2, 2> m) {
  auto roots_obj = get_eigenvalues(m);
  auto nroots = roots_obj.n;
  auto roots = roots_obj.values;
  auto mults = roots_obj.mults;
  /* there are only a few output cases, see solve_quadratic() */
  Matrix<2, 2> q;
  Vector<2> l;
  if (nroots == 2) {
    for (Int i = 0; i < 2; ++i) {
      q[i] = single_eigenvector(m, roots[i]);
      l[i] = roots[i];
    }
  } else {
    CHECK(nroots == 1 && mults[0] == 2);
    l[0] = l[1] = roots[0];
    q = identity_matrix<2, 2>();
  }
  return {q, l};
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
INLINE DiagDecomp<dim> decompose_eigen(Matrix<dim, dim> m) {
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
   are the right eigenvectors, *not* the
   change of basis matrix */
template <Int dim>
INLINE Matrix<dim, dim> compose_eigen(Matrix<dim, dim> q, Vector<dim> l) {
  return transpose(q * diagonal(l) * invert(q));
}

/* like the above, but knowing Q is orthonormal,
   meaning in this case it *is* the change of basis
   matrix */
template <Int dim>
INLINE Matrix<dim, dim> compose_ortho(Matrix<dim, dim> q, Vector<dim> l) {
  return q * diagonal(l) * transpose(q);
}
Reals get_max_eigenvalues(Int dim, Reals symms);

}  // end namespace Omega_h

#endif
