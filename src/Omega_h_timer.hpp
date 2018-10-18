#ifndef OMEGA_H_TIMER_HPP
#define OMEGA_H_TIMER_HPP

#include <chrono>

#include <Omega_h_defines.hpp>

namespace Omega_h {

struct Now {
  typedef std::chrono::time_point<std::chrono::steady_clock> Impl;
  Impl impl;
};

Now now();

Real operator-(Now b, Now a);

}  // end namespace Omega_h

#endif
