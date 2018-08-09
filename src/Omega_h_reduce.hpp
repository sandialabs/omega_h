#ifndef OMEGA_H_REDUCE_HPP
#define OMEGA_H_REDUCE_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_stack.hpp>

#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_kokkos.hpp>
#endif

namespace Omega_h {

template <typename T>
typename T::value_type parallel_reduce(LO n, T f, char const* name = "") {
  using VT = typename T::value_type;
#ifdef OMEGA_H_USE_KOKKOSCORE
  VT result;
  f.init(result);
  if (n > 0) Kokkos::parallel_reduce(name, policy(n), f, result);
#else
  begin_code(name);
  VT result;
  f.init(result);
  for (LO i = 0; i < n; ++i) f(i, result);
  end_code();
#endif
  return result;
}

}

#endif
