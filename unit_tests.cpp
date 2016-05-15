#include "internal.hpp"

template <UInt m, UInt n>
static void test_qr_decomp(Matrix<m,n> a) {
  Matrix<m,n> q;
  Matrix<n,n> r;
  decompose_qr_reduced(a, q, r);
  CHECK(are_close(a, q * r));
  CHECK(are_close(transpose(q) * q, identity_matrix<m,n>()));
}

static void test_qr_decomps() {
  test_qr_decomp(Matrix<3,3>({
       0, 0, 0,
       0, 0, 0,
       0, 0, 0}));
  test_qr_decomp(Matrix<3,3>({
       EPSILON, 0, 0,
       EPSILON, EPSILON, 0,
       EPSILON, EPSILON, EPSILON}));
  test_qr_decomp(Matrix<3,3>({
      12, -51,  4,
       6, 167,-68,
      -4,  24,-41}));
}

static void test_form_ortho_basis() {
  auto n = normalize(vector_3(1,1,1));
  auto f = form_ortho_basis(n);
  CHECK(are_close(f[0], n));
  CHECK(are_close(transpose(f) * f, identity_matrix<3,3>()));
}

static void test_least_squares() {
  Matrix<4,2> m({
      1, 1,
      1, 2,
      1, 3,
      1, 4});
  Vector<4> b({6, 5, 7, 10});
  Vector<2> x;
  CHECK(solve_least_squares_qr(m, b, x));
  CHECK(are_close(x, vector_2(3.5, 1.4)));
}

static void test_int128() {
  Int128 a(INT64_MAX);
  auto b = a + a;
  b = b + b;
  b = b + b;
  b = b >> 3;
  CHECK(b == a);
}

static void test_repro_sum() {
  Reals a({std::exp2(int(20)),std::exp2(int(-20))});
  Real sum = repro_sum(a);
  CHECK(sum == std::exp2(20) + std::exp2(int(-20)));
}

static void test_cubic(Real a, Real b, Real c,
    UInt nroots_wanted, Few<Real, 3> roots_wanted,
    Few<UInt, 3> mults_wanted) {
  Few<Real, 3> roots;
  Few<UInt, 3> mults;
  UInt nroots = solve_cubic(a, b, c, roots, mults);
  CHECK(nroots == nroots_wanted);
  for (UInt i = 0; i < nroots; ++i) {
    CHECK(mults[i] == mults_wanted[i]);
    CHECK(are_close(roots[i], roots_wanted[i]));
  }
}

static void test_cubic() {
  test_cubic(0, 0, 0,
      1, Few<Real,3>({0}), Few<UInt,3>({3}));
  test_cubic(-3. / 2., -3. / 2., 1.,
      3, Few<Real,3>({-1,2,.5}), Few<UInt,3>({1,1,1}));
  test_cubic(0, -3., 2.,
      2, Few<Real,3>({-2,1}), Few<UInt,3>({1,2}));
  test_cubic(3, -6, -8,
      3, Few<Real,3>({-4,2,-1}), Few<UInt,3>({1,1,1}));
}

static void test_eigen_cubic(Matrix<3,3> m,
    Matrix<3,3> q_expect, Vector<3> l_expect) {
  Matrix<3,3> q;
  Vector<3> l;
  bool ok = decompose_eigen_polynomial(m, q, l);
  CHECK(ok);
  CHECK(are_close(q,q_expect));
  CHECK(are_close(l,l_expect));
}

static void test_eigen_cubic(Matrix<3,3> m,
    Vector<3> l_expect) {
  Matrix<3,3> q;
  Vector<3> l;
  bool ok = decompose_eigen_polynomial(m, q, l);
  CHECK(ok);
  CHECK(are_close(l,l_expect,1e-6));
  CHECK(are_close(m, transpose(q * diagonal(l) * invert(q))));
}

static void test_eigen_cubic_ortho(Matrix<3,3> m,
    Vector<3> l_expect) {
  Matrix<3,3> q;
  Vector<3> l;
  bool ok = decompose_eigen_polynomial(m, q, l);
  CHECK(ok);
  CHECK(are_close(transpose(q) * q, identity_matrix<3,3>()));
  CHECK(are_close(l,l_expect));
  CHECK(are_close(m, (q * diagonal(l) * transpose(q))));
}

static void test_eigen_cubic() {
  test_eigen_cubic(
      identity_matrix<3,3>(),
      identity_matrix<3,3>(),
      vector_3(1,1,1));
  test_eigen_cubic(
      matrix_3x3(0,0,0,0,0,0,0,0,0),
      identity_matrix<3,3>(),
      vector_3(0,0,0));
  test_eigen_cubic(
      matrix_3x3(
        -1, 3, -1,
        -3, 5, -1,
        -3, 3,  1),
      vector_3(1,2,2));
  {
  auto q = rotate(PI / 4, vector_3(0,0,1)) *
           rotate(PI / 4, vector_3(0,1,0));
  CHECK(are_close(transpose(q) * q, identity_matrix<3,3>()));
  auto l = diagonal(vector_3(1, 1, 1e-6));
  auto a = q * l * transpose(q);
  test_eigen_cubic_ortho(a, vector_3(1e-6, 1, 1));
  }
  {
  auto q = rotate(PI / 4, vector_3(0,0,1)) *
           rotate(PI / 4, vector_3(0,1,0));
  CHECK(are_close(transpose(q) * q, identity_matrix<3,3>()));
  auto l = diagonal(vector_3(1e-6, 1e-6, 1));
  auto a = q * l * transpose(q);
  test_eigen_cubic_ortho(a, vector_3(1, 1e-6, 1e-6));
  }
}

int main(int argc, char** argv) {
  init(argc, argv);
  test_form_ortho_basis();
  test_qr_decomps();
  test_eigen_cubic();
  test_least_squares();
  test_int128();
  test_repro_sum();
  test_cubic();
  fini();
}
