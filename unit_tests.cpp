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

  auto q = rotate(PI / 4, vector_3(0,0,1)) *
           rotate(PI / 4, vector_3(0,1,0));
  std::cout << "q\n" << q;
  CHECK(are_close(transpose(q) * q, identity_matrix<3,3>()));
  auto l = matrix_3x3(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1e-6);
  std::cout << "l\n" << l;
  auto a = q * l * transpose(q);
  std::cout << "a\n" << a;
  Matrix<3,3> q2;
  Matrix<3,3> l2;
  decompose_eigen_qr(a, q2, l2);
  std::cout << "q2\n" << q2;
  std::cout << "l2\n" << l2;
}
