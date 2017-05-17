#include "Omega_h_int128.hpp"

namespace Omega_h {

double Int128::to_double(double unit) const {
  Int128 tmp = *this;
  if (tmp < Int128(0)) tmp = -tmp;
  while (tmp.high) {
    tmp = tmp >> 1;
    unit *= 2;
  }
  double x = tmp.low;
  if (*this < Int128(0)) x = -x;
  x *= unit;
  return x;
}

}  // end namespace Omega_h
