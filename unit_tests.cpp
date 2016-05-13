#include "internal.hpp"

template <UInt m, UInt n>
static void test_qr_decomp(Matrix<m,n> a) {
  Matrix<m,n> q;
  Matrix<n,n> r;
  decompose_qr_reduced(a, q, r);
  CHECK(are_close(a, q * r));
  CHECK(are_close(transpose(q) * q, identity_matrix<m,n>()));
}

int main() {
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

  Matrix<3,3> a({
      1, 0, 0,
      0, 1, 0,
      0, 0, 1});
  Matrix<3,3> q;
  Matrix<3,3> l;
  decompose_eigen_qr(a, q, l, 3);
  std::cout << q << '\n';
  std::cout << l << '\n';
}
