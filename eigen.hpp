#ifndef EIGEN_HPP
#define EIGEN_HPP

namespace osh {

/* http://mathworld.wolfram.com/CharacteristicPolynomial.html */
INLINE void characteristic_cubic(Matrix<3, 3> A, Real& a, Real& b, Real& c) {
  Real tA = trace(A);
  Real c2 = (1. / 2.) * ((tA * tA) - trace(A * A));
  a = -tA;
  b = c2;
  c = -determinant(A);
}

/* http://mathworld.wolfram.com/CharacteristicPolynomial.html */
INLINE void characteristic_quadratic(Matrix<2, 2> A, Real& a, Real& b) {
  a = -trace(A);
  b = determinant(A);
}

/* the null space of the matrix (s = m - l*I)
   is the space spanned by the eigenvectors.
   the multiplicity of this root is the dimensionality
   of that space, i.e. the number of eigenvectors. */

/* in the case that the null space is 1D and space is 3D,
   take the largest cross product of any pair of columns */
INLINE void single_eigenvector(Matrix<3, 3> m, Real l, Vector<3>& v) {
  subtract_from_diag(m, l);
  v = cross(m[0], m[1]);
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
INLINE void double_eigenvector(Matrix<3, 3> m, Real l, Vector<3>& u,
                               Vector<3>& v) {
  subtract_from_diag(m, l);
  Vector<3> n = get_1d_column_space(m);
  Matrix<3, 3> b = form_ortho_basis(n);
  u = b[1];
  v = b[2];
}

INLINE void decompose_eigen2(Matrix<3, 3> m, Matrix<3, 3>& q, Vector<3>& l) {
  Real a, b, c;
  characteristic_cubic(m, a, b, c);
  Few<Real, 3> roots;
  Few<Int, 3> mults;
  Int nroots = solve_cubic(a, b, c, roots, mults);
  /* there are only a few output cases, see solve_cubic() */
  if (nroots == 3) {
    for (Int i = 0; i < 3; ++i) {
      single_eigenvector(m, roots[i], q[i]);
      l[i] = roots[i];
    }
  } else if (nroots == 2 && mults[1] == 2) {
    single_eigenvector(m, roots[0], q[0]);
    l[0] = roots[0];
    double_eigenvector(m, roots[1], q[1], q[2]);
    l[1] = l[2] = roots[1];
  } else {
    CHECK(nroots == 1 && mults[0] == 3);
    l[0] = l[1] = l[2] = roots[0];
    q = identity_matrix<3, 3>();
  }
}

/* in the case that the null space is 1D and space is 2D,
   get the largest column and rotate it 90 deg */
INLINE void single_eigenvector(Matrix<2, 2> m, Real l, Vector<2>& v) {
  Matrix<2, 2> s = (m - (l * identity_matrix<2, 2>()));
  v = perp(get_1d_column_space(s));
}

INLINE void decompose_eigen2(Matrix<2, 2> m, Matrix<2, 2>& q, Vector<2>& l) {
  Real a, b;
  characteristic_quadratic(m, a, b);
  Few<Real, 2> roots;
  Few<Int, 2> mults;
  Int nroots = solve_quadratic(a, b, roots, mults);
  /* there are only a few output cases, see solve_quadratic() */
  if (nroots == 2) {
    for (Int i = 0; i < 2; ++i) {
      single_eigenvector(m, roots[i], q[i]);
      l[i] = roots[i];
    }
  } else {
    CHECK(nroots == 1 && mults[0] == 2);
    l[0] = l[1] = roots[0];
    q = identity_matrix<2, 2>();
  }
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
INLINE void decompose_eigen(Matrix<dim, dim> m, Matrix<dim, dim>& q,
                            Vector<dim>& l) {
  /* the cubic solver is especially sensitive to dynamic
     range. what we can do is to normalize the input matrix
     and then re-apply that norm to the resulting roots */
  Real nm = max_norm(m);
  if (nm > EPSILON) {
    m = m / nm;
    decompose_eigen2(m, q, l);
    l = l * nm;
  } else {
    /* this is the zero matrix... */
    q = identity_matrix<dim, dim>();
    l = zero_vector<dim>();
  }
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

} //end namespace osh

#endif
