#ifndef OMEGA_H_MATRIX_HPP
#define OMEGA_H_MATRIX_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_vector.hpp>

namespace Omega_h {

/* column-first storage and indexing !!! */

#ifdef OMEGA_H_USE_KOKKOS

template <Int m, Int n>
class Matrix : public Few<Vector<m>, n> {
 public:
  OMEGA_H_INLINE Matrix() {}
  /* these constructors accept the matrix in
   * row-first order for convenience.
   */
  inline Matrix(std::initializer_list<Vector<m>> l) : Few<Vector<m>, n>(l) {}
  inline Matrix(std::initializer_list<Real> l);
  OMEGA_H_INLINE void operator=(Matrix<m, n> const& rhs) volatile {
    Few<Vector<m>, n>::operator=(rhs);
  }
  OMEGA_H_INLINE Matrix(Matrix<m, n> const& rhs) : Few<Vector<m>, n>(rhs) {}
  OMEGA_H_INLINE Matrix(const volatile Matrix<m, n>& rhs)
      : Few<Vector<m>, n>(rhs) {}
  OMEGA_H_INLINE Real& operator()(Int i, Int j) {
    return Few<Vector<m>, n>::operator[](j)[i];
  }
  OMEGA_H_INLINE Real const& operator()(Int i, Int j) const {
    return Few<Vector<m>, n>::operator[](j)[i];
  }
  OMEGA_H_INLINE Real volatile& operator()(Int i, Int j) volatile {
    return Few<Vector<m>, n>::operator[](j)[i];
  }
  OMEGA_H_INLINE Real const volatile& operator()(Int i, Int j) const volatile {
    return Few<Vector<m>, n>::operator[](j)[i];
  }
};

#else

template <Int m, Int n>
class Matrix : public Few<Vector<m>, n> {
 public:
  inline Matrix() = default;
  /* these constructors accept the matrix in
   * row-first order for convenience.
   */
  inline Matrix(std::initializer_list<Vector<m>> l) : Few<Vector<m>, n>(l) {}
  inline Matrix(std::initializer_list<Real> l);
  inline Matrix(Matrix const&) = default;
  inline Matrix(Matrix&&) = default;
  inline Matrix& operator=(Matrix const&) = default;
  inline Matrix& operator=(Matrix&&) = default;
  OMEGA_H_INLINE Real& operator()(Int i, Int j) OMEGA_H_NOEXCEPT {
    return Few<Vector<m>, n>::operator[](j)[i];
  }
  OMEGA_H_INLINE Real const& operator()(Int i, Int j) const OMEGA_H_NOEXCEPT {
    return Few<Vector<m>, n>::operator[](j)[i];
  }
};

#endif

template <Int dim>
using Tensor = Matrix<dim, dim>;

template <Int m, Int n>
OMEGA_H_INLINE Real* scalar_ptr(Matrix<m, n>& a) {
  return &a[0][0];
}
template <Int m, Int n>
OMEGA_H_INLINE Real const* scalar_ptr(Matrix<m, n> const& a) {
  return &a[0][0];
}

template <Int m, Int n>
inline Matrix<m, n>::Matrix(std::initializer_list<Real> l) {
  Int k = 0;
  for (Real v : l) {
    (*this)[k % n][k / n] = v;
    ++k;
  }
}

template <Int m, Int n>
OMEGA_H_INLINE Vector<m> operator*(
    Matrix<m, n> a, Vector<n> b)OMEGA_H_NOEXCEPT {
  Vector<m> c = a[0] * b[0];
  for (Int j = 1; j < n; ++j) c = c + a[j] * b[j];
  return c;
}

template <Int m, Int n, Int p>
OMEGA_H_INLINE Matrix<m, n> operator*(Matrix<m, p> a, Matrix<p, n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a * b[j];
  return c;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<n, m> transpose(Matrix<m, n> a) {
  Matrix<n, m> b;
  for (Int i = 0; i < m; ++i)
    for (Int j = 0; j < n; ++j) b[i][j] = a[j][i];
  return b;
}

template <Int max_m, Int max_n>
OMEGA_H_INLINE Matrix<max_m, max_n> identity_matrix(Int m, Int n) {
  Matrix<max_m, max_n> a;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < m; ++i) a[j][i] = (i == j);
  return a;
}

template <Int max_m, Int max_n>
OMEGA_H_INLINE Matrix<max_m, max_n> identity_matrix() {
  return identity_matrix<max_m, max_n>(max_m, max_n);
}

template <Int dim>
OMEGA_H_INLINE Tensor<dim> identity_tensor() {
  return identity_matrix<dim, dim>();
}

OMEGA_H_INLINE Matrix<1, 1> matrix_1x1(Real a) {
  Matrix<1, 1> o;
  o[0][0] = a;
  return o;
}

OMEGA_H_INLINE Tensor<1> tensor_1(Real a) { return matrix_1x1(a); }

OMEGA_H_INLINE Matrix<2, 2> matrix_2x2(Real a, Real b, Real c, Real d) {
  Matrix<2, 2> o;
  o[0] = vector_2(a, c);
  o[1] = vector_2(b, d);
  return o;
}

OMEGA_H_INLINE Tensor<2> tensor_2(Real a, Real b, Real c, Real d) {
  return matrix_2x2(a, b, c, d);
}

OMEGA_H_INLINE Matrix<3, 3> matrix_3x3(
    Real a, Real b, Real c, Real d, Real e, Real f, Real g, Real h, Real i) {
  Matrix<3, 3> o;
  o[0] = vector_3(a, d, g);
  o[1] = vector_3(b, e, h);
  o[2] = vector_3(c, f, i);
  return o;
}

OMEGA_H_INLINE Tensor<3> tensor_3(
    Real a, Real b, Real c, Real d, Real e, Real f, Real g, Real h, Real i) {
  return matrix_3x3(a, b, c, d, e, f, g, h, i);
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> operator*(Matrix<m, n> a, Real b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a[j] * b;
  return c;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n>& operator*=(Matrix<m, n>& a, Real b) {
  a = a * b;
  return a;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> operator*(Real a, Matrix<m, n> b) {
  return b * a;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> operator/(Matrix<m, n> a, Real b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a[j] / b;
  return c;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n>& operator/=(Matrix<m, n>& a, Real b) {
  a = a / b;
  return a;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> operator+(Matrix<m, n> a, Matrix<m, n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a[j] + b[j];
  return c;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n>& operator+=(Matrix<m, n>& a, Matrix<m, n> b) {
  a = a + b;
  return a;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> operator-(Matrix<m, n> a, Matrix<m, n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a[j] - b[j];
  return c;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n>& operator-=(Matrix<m, n>& a, Matrix<m, n> b) {
  a = a - b;
  return a;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> operator-(Matrix<m, n> a) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = -a[j];
  return c;
}

template <Int m, Int n>
OMEGA_H_INLINE Real max_norm(Matrix<m, n> a) {
  Real x = 0.0;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < m; ++i) x = max2(x, std::abs(a[j][i]));
  return x;
}

template <Int m, Int n>
OMEGA_H_INLINE bool are_close(
    Matrix<m, n> a, Matrix<m, n> b, Real tol = EPSILON, Real floor = EPSILON) {
  for (Int j = 0; j < n; ++j)
    if (!are_close(a[j], b[j], tol, floor)) return false;
  return true;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> outer_product(Vector<m> a, Vector<n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < m; ++i) c[j][i] = a[i] * b[j];
  return c;
}

template <Int m>
OMEGA_H_INLINE Real trace(Tensor<m> a) {
  Real t = a[0][0];
  for (Int i = 1; i < m; ++i) t += a[i][i];
  return t;
}

template <Int m>
OMEGA_H_INLINE Vector<m> diagonal(Tensor<m> a) {
  Vector<m> v;
  for (Int i = 0; i < m; ++i) v[i] = a[i][i];
  return v;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> zero_matrix() {
  Matrix<m, n> a;
  for (Int j = 0; j < n; ++j) a[j] = zero_vector<m>();
  return a;
}

OMEGA_H_INLINE Real determinant(Tensor<1> m) { return m[0][0]; }

OMEGA_H_INLINE Real determinant(Tensor<2> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[0][1];
  Real d = m[1][1];
  return a * d - b * c;
}

OMEGA_H_INLINE Real determinant(Tensor<3> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[2][0];
  Real d = m[0][1];
  Real e = m[1][1];
  Real f = m[2][1];
  Real g = m[0][2];
  Real h = m[1][2];
  Real i = m[2][2];
  return (a * e * i) + (b * f * g) + (c * d * h) - (c * e * g) - (b * d * i) -
         (a * f * h);
}

template <Int max_m, Int max_n>
OMEGA_H_INLINE Real inner_product(
    Int m, Int n, Matrix<max_m, max_n> a, Matrix<max_m, max_n> b) {
  Real out = 0.0;
  for (Int j = 0; j < n; ++j) {
    for (Int i = 0; i < m; ++i) {
      out += a[j][i] * b[j][i];
    }
  }
  return out;
}

template <Int m, Int n>
OMEGA_H_INLINE Real inner_product(Matrix<m, n> a, Matrix<m, n> b) {
  return inner_product(m, n, a, b);
}

template <Int max_m, Int max_n>
OMEGA_H_INLINE Real norm(
    Int const m, Int const n, Matrix<max_m, max_n> const a) {
  return std::sqrt(inner_product(m, n, a, a));
}

template <Int m, Int n>
OMEGA_H_INLINE Real norm(Matrix<m, n> const a) {
  return norm(m, n, a);
}

OMEGA_H_INLINE Tensor<3> cross(Vector<3> const a) {
  return matrix_3x3(0, -a[2], a[1], a[2], 0, -a[0], -a[1], a[0], 0);
}

OMEGA_H_INLINE Vector<3> uncross(Tensor<3> const c) {
  return 0.5 *
         vector_3(c[1][2] - c[2][1], c[2][0] - c[0][2], c[0][1] - c[1][0]);
}

OMEGA_H_INLINE Tensor<1> invert(Tensor<1> const m) {
  return matrix_1x1(1.0 / m[0][0]);
}

OMEGA_H_INLINE Tensor<2> invert(Tensor<2> const m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[0][1];
  Real d = m[1][1];
  return matrix_2x2(d, -b, -c, a) / determinant(m);
}

OMEGA_H_INLINE Tensor<3> invert(Tensor<3> const a) {
  Tensor<3> b;
  b[0] = cross(a[1], a[2]);
  b[1] = cross(a[2], a[0]);
  b[2] = cross(a[0], a[1]);
  return transpose(b) / determinant(a);
}

template <Int m, Int n>
OMEGA_H_INLINE typename std::enable_if<(n < m), Matrix<n, m>>::type
pseudo_invert(Matrix<m, n> a) {
  auto at = transpose(a);
  return invert(at * a) * at;
}

template <Int m, Int n>
OMEGA_H_INLINE typename std::enable_if<(n > m), Matrix<n, m>>::type
pseudo_invert(Matrix<m, n> a) {
  auto at = transpose(a);
  return at * invert(a * at);
}

template <Int m>
OMEGA_H_INLINE Tensor<m> diagonal(Vector<m> v) {
  Tensor<m> a;
  for (Int i = 0; i < m; ++i)
    for (Int j = i + 1; j < m; ++j) a[i][j] = a[j][i] = 0.0;
  for (Int i = 0; i < m; ++i) a[i][i] = v[i];
  return a;
}

OMEGA_H_INLINE constexpr Int symm_ncomps(Int dim) {
  return ((dim + 1) * dim) / 2;
}

OMEGA_H_INLINE Vector<1> symm2vector(Tensor<1> symm) {
  return vector_1(symm[0][0]);
}

OMEGA_H_INLINE Tensor<1> vector2symm(Vector<1> v) { return matrix_1x1(v[0]); }

OMEGA_H_INLINE Vector<3> symm2vector(Tensor<2> symm) {
  Vector<3> v;
  v[0] = symm[0][0];
  v[1] = symm[1][1];
  v[2] = symm[1][0];
  return v;
}

OMEGA_H_INLINE Tensor<2> vector2symm(Vector<3> v) {
  Tensor<2> symm;
  symm[0][0] = v[0];
  symm[1][1] = v[1];
  symm[1][0] = v[2];
  symm[0][1] = symm[1][0];
  return symm;
}

OMEGA_H_INLINE Vector<6> symm2vector(Tensor<3> symm) {
  Vector<6> v;
  v[0] = symm[0][0];
  v[1] = symm[1][1];
  v[2] = symm[2][2];
  v[3] = symm[1][0];
  v[4] = symm[2][1];
  v[5] = symm[2][0];
  return v;
}

OMEGA_H_INLINE Tensor<3> vector2symm(Vector<6> v) {
  Tensor<3> symm;
  symm[0][0] = v[0];
  symm[1][1] = v[1];
  symm[2][2] = v[2];
  symm[1][0] = v[3];
  symm[2][1] = v[4];
  symm[2][0] = v[5];
  symm[0][1] = symm[1][0];
  symm[1][2] = symm[2][1];
  symm[0][2] = symm[2][0];
  return symm;
}

/* Symmetric metric tensor storage convention used by
   INRIA:  https://hal.inria.fr/inria-00363007/document */
OMEGA_H_INLINE Tensor<1> vector2symm_inria(Vector<1> v) {
  return matrix_1x1(v[0]);
}

OMEGA_H_INLINE Tensor<2> vector2symm_inria(Vector<3> v) {
  Tensor<2> symm;
  symm[0][0] = v[0];
  symm[0][1] = v[1];
  symm[1][1] = v[2];
  symm[1][0] = symm[0][1];
  return symm;
}

OMEGA_H_INLINE Tensor<3> vector2symm_inria(Vector<6> v) {
  Tensor<3> symm;
  symm[0][0] = v[0];
  symm[0][1] = v[1];
  symm[1][1] = v[2];
  symm[0][2] = v[3];
  symm[1][2] = v[4];
  symm[2][2] = v[5];
  symm[1][0] = symm[0][1];
  symm[2][0] = symm[0][2];
  symm[2][1] = symm[1][2];
  return symm;
}

OMEGA_H_INLINE Vector<1> symm2vector_inria(Tensor<1> symm) {
  return vector_1(symm[0][0]);
}

OMEGA_H_INLINE Vector<3> symm2vector_inria(Tensor<2> symm) {
  Vector<3> v;
  v[0] = symm[0][0];
  v[1] = symm[0][1];
  v[2] = symm[1][1];
  return v;
}

OMEGA_H_INLINE Vector<6> symm2vector_inria(Tensor<3> symm) {
  Vector<6> v;
  v[0] = symm[0][0];
  v[1] = symm[0][1];
  v[2] = symm[1][1];
  v[3] = symm[0][2];
  v[4] = symm[1][2];
  v[5] = symm[2][2];
  return v;
}

OMEGA_H_INLINE constexpr Int matrix_ncomps(Int m, Int n) { return m * n; }

template <Int m, Int n>
OMEGA_H_INLINE Vector<matrix_ncomps(m, n)> matrix2vector(Matrix<m, n> a) {
  Vector<matrix_ncomps(m, n)> v;
  for (Int i = 0; i < m; ++i) {
    for (Int j = 0; j < n; ++j) {
      v[i * n + j] = a(i, j);
    }
  }
  return v;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<m, n> vector2matrix(Vector<matrix_ncomps(m, n)> v) {
  Matrix<m, n> a;
  for (Int i = 0; i < m; ++i) {
    for (Int j = 0; j < n; ++j) {
      a(i, j) = v[i * n + j];
    }
  }
  return a;
}

template <Int n>
OMEGA_H_DEVICE void set_symm(Write<Real> const& a, Int i, Tensor<n> symm) {
  set_vector(a, i, symm2vector(symm));
}

template <Int n, typename Arr>
OMEGA_H_DEVICE Tensor<n> get_symm(Arr const& a, Int i) {
  return vector2symm(get_vector<symm_ncomps(n)>(a, i));
}

template <Int m, Int n>
OMEGA_H_DEVICE void set_matrix(
    Write<Real> const& matrices, Int i, Matrix<m, n> matrix) {
  set_vector(matrices, i, matrix2vector(matrix));
}

template <Int m, Int n, typename Arr>
OMEGA_H_DEVICE Matrix<m, n> get_matrix(Arr const& matrices, Int const i) {
  return vector2matrix<m, n>(get_vector<matrix_ncomps(m, n)>(matrices, i));
}

/* Rodrigues' Rotation Formula */
OMEGA_H_INLINE Tensor<3> rotate(Real const angle, Vector<3> const axis) {
  return std::cos(angle) * identity_matrix<3, 3>() +
         std::sin(angle) * cross(axis) +
         (1 - std::cos(angle)) * outer_product(axis, axis);
}

OMEGA_H_INLINE Tensor<2> rotate(Real const angle) {
  return matrix_2x2(
      std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle));
}

OMEGA_H_INLINE Real rotation_angle(Tensor<2> const r) {
  auto const cos_theta = 0.5 * trace(r);
  auto const sin_theta = 0.5 * (r(1, 0) - r(0, 1));
  return std::atan2(sin_theta, cos_theta);
}

OMEGA_H_INLINE Real rotation_angle(Tensor<3> const r) {
  auto const cos_theta = 0.5 * (trace(r) - 1.0);
  return std::acos(cos_theta);
}

OMEGA_H_INLINE Tensor<1> form_ortho_basis(Vector<1> const v) {
  return tensor_1(v[0]);
}

OMEGA_H_INLINE Tensor<2> form_ortho_basis(Vector<2> const v) {
  Tensor<2> A;
  A[0] = v;
  A[1] = perp(A[0]);
  return A;
}

/* Duff, Tom, et al.
   "Building an orthonormal basis, revisited."
   Journal of Computer Graphics Techniques Vol 6.1 (2017). */
OMEGA_H_INLINE Tensor<3> form_ortho_basis(Vector<3> const v) {
  Tensor<3> A;
  A[0] = v;
  auto sign = std::copysign(1.0, v(2));
  const auto a = -1.0 / (sign + v(2));
  const auto b = v(0) * v(1) * a;
  A[1] = vector_3(1.0 + sign * v(0) * v(0) * a, sign * b, -sign * v(0));
  A[2] = vector_3(b, sign + v(1) * v(1) * a, -v(1));
  return A;
}

template <Int dim>
OMEGA_H_INLINE Tensor<dim> deviator(Tensor<dim> const a) {
  return (a - ((1.0 / dim) * trace(a) * identity_matrix<dim, dim>()));
}

template <int new_dim, int old_dim>
OMEGA_H_INLINE Tensor<new_dim> resize(Tensor<old_dim> const m) {
  constexpr int min_dim = Omega_h::min2(new_dim, old_dim);
  Tensor<new_dim> m2;
  for (int i = 0; i < min_dim; ++i)
    for (int j = 0; j < min_dim; ++j) m2(i, j) = m(i, j);
  for (int i = 0; i < new_dim; ++i)
    for (int j = min_dim; j < new_dim; ++j) m2(i, j) = m2(j, i) = 0.0;
  return m2;
}

template <Int dim>
Reals repeat_symm(LO const n, Tensor<dim> const symm);
extern template Reals repeat_symm(LO const n, Tensor<3> const symm);
extern template Reals repeat_symm(LO const n, Tensor<2> const symm);
extern template Reals repeat_symm(LO const n, Tensor<1> const symm);

Reals resize_symms(Reals old_symms, Int old_dim, Int new_dim);

template <Int dim>
Reals repeat_matrix(LO const n, Tensor<dim> m);

extern template Reals repeat_matrix(LO const n, Tensor<3> const m);
extern template Reals repeat_matrix(LO const n, Tensor<2> const m);
extern template Reals repeat_matrix(LO const n, Tensor<1> const m);

Reals matrices_times_vectors(Reals ms, Reals vs, Int dim);
Reals matrices_times_matrices(Reals ms, Reals vs, Int dim);

Reals symms_inria2osh(Int dim, Reals symms);
Reals symms_osh2inria(Int dim, Reals symms);

Reals matrices_to_symms(Reals const matrices, Int const dim);

}  // end namespace Omega_h

#endif
