#include "internal.hpp"

int main()
{
  Matrix<3,3> a({
      12, -51,  4,
       6, 167,-68,
      -4,  24,-41});
//     1,   0,  0,
//     0,   1,  0,
//     0,   0,  1});
  Matrix<3,3> q;
  Matrix<3,3> r;
  decompose_qr_reduced(a, q, r);
  std::cout << a << '\n';
  std::cout << q << '\n';
  std::cout << r << '\n';
  std::cout << (q * r) << '\n';
  CHECK(are_close(a, q * r, 1e-12, 1e-10));
}
