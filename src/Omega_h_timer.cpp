#include "Omega_h_timer.hpp"

namespace Omega_h {

Now now() {
  Now t;
  t.impl = std::chrono::steady_clock::now();
  return t;
}

Real operator-(Now b, Now a) {
  return std::chrono::duration<double>(b.impl - a.impl).count();
}

}  // end namespace Omega_h
