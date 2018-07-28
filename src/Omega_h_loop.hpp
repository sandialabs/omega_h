#ifndef OMEGA_H_LOOP_HPP
#define OMEGA_H_LOOP_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_stack.hpp>

namespace Omega_h {

#ifdef OMEGA_H_USE_KOKKOSCORE
using ExecSpace = Kokkos::DefaultExecutionSpace;
using StaticSched = Kokkos::Schedule<Kokkos::Static>;
using Policy = Kokkos::RangePolicy<ExecSpace, StaticSched, Omega_h::LO>;

inline Policy policy(LO n) { return Policy(0, n); }
#endif

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
