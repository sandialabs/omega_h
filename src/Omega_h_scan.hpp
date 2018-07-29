#ifndef OMEGA_H_LOOP_HPP
#define OMEGA_H_LOOP_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_stack.hpp>

#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_kokkos.hpp>
#endif

namespace Omega_h {

template <typename T>
void parallel_scan(LO n, T f, char const* name = "") {
#ifdef OMEGA_H_USE_KOKKOSCORE
  if (n > 0) Kokkos::parallel_scan(name, policy(n), f);
#else
  using VT = typename T::value_type;
  begin_code(name);
  VT update;
  f.init(update);
  for (LO i = 0; i < n; ++i) f(i, update, true);
  end_code();
#endif
}

}  // end namespace Omega_h

#endif

