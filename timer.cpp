#include "timer.hpp"

namespace osh {

Now now() {
  Now t;
  t.impl = std::chrono::high_resolution_clock::now();
  return t;
}

Real operator-(Now b, Now a) {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(b.impl - a.impl)
             .count() *
         1e-9;
}

} //end namespace osh
