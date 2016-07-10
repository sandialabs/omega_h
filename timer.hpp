#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

#include "internal.hpp"

namespace osh {

struct Now {
  typedef std::chrono::time_point<std::chrono::high_resolution_clock> Impl;
  Impl impl;
};

Now now();

Real operator-(Now b, Now a);

}  // end namespace osh

#endif
