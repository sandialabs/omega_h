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

static void test_eigen_decomp() {
  auto q = rotate(PI / 4, vector_3(0,0,1)) *
           rotate(PI / 4, vector_3(0,1,0));
  CHECK(are_close(transpose(q) * q, identity_matrix<3,3>()));
  auto l = matrix_3x3(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1e-6);
  auto a = q * l * transpose(q);
  Matrix<3,3> q2;
  Matrix<3,3> l2;
  decompose_eigen_qr(a, q2, l2);
  CHECK(are_close(transpose(q2) * q2, identity_matrix<3,3>()));
  CHECK(are_close(q2 * l2 * transpose(q2), a));
}

static void test_least_squares() {
  Matrix<4,2> m({
      1, 1,
      1, 2,
      1, 3,
      1, 4});
  Vector<4> b({6, 5, 7, 10});
  Vector<2> x = solve_least_squares_qr(m, b);
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

int main(int argc, char** argv) {
  init(argc, argv);
  test_qr_decomps();
  test_eigen_decomp();
  test_least_squares();
  test_int128();
  test_repro_sum();
  fini();
}
