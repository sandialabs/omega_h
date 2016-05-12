#include "internal.hpp"

int main()
{
  Matrix<3,3> a({
      12, -51,  4,
       6, 167,-68,
      -4,  24,-41});
  Matrix<3,3> q;
  Matrix<3,3> r;
  decompose_qr_reduced(a, q, r);
}
