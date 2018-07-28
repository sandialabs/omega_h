#ifndef OMEGA_H_FOR_HPP
#define OMEGA_H_FOR_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_stack.hpp>

#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_kokkos.hpp>
#endif

namespace Omega_h {

template <typename T>
void parallel_for(LO n, T const& f, char const* name = "") {
#ifdef OMEGA_H_USE_KOKKOSCORE
  if (n > 0) Kokkos::parallel_for(name, policy(n), f);
#else
  begin_code(name);
  for (LO i = 0; i < n; ++i) f(i);
  end_code();
#endif
}

}  // end namespace Omega_h

#endif
