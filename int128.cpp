#include "int128.hpp"

namespace osh {

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

void Int128::print(std::ostream& o) const {
  std::ios::fmtflags f(o.flags());
  o << std::hex << "int128(" << high << ',' << low << ')';
  o.flags(f);
}

std::ostream& operator<<(std::ostream& o, Int128 const& x) {
  x.print(o);
  return o;
}

} //end namespace osh
