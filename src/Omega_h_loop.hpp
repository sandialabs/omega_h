#ifndef OMEGA_H_LOOP_HPP
#define OMEGA_H_LOOP_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_kokkos.hpp>

namespace Omega_h {

#ifdef OMEGA_H_USE_KOKKOSCORE
using ExecSpace = Kokkos::DefaultExecutionSpace;
using StaticSched = Kokkos::Schedule<Kokkos::Static>;
using Policy = Kokkos::RangePolicy<ExecSpace, StaticSched>;

inline Policy policy(LO n) { return Policy(0, static_cast<std::size_t>(n)); }
#endif

template <typename T>
void parallel_for(LO n, T const& f, std::string const& name = "") {
#ifdef OMEGA_H_USE_KOKKOSCORE
  if (n > 0) Kokkos::parallel_for(policy(n), f, name);
#else
  begin_code(name);
  for (LO i = 0; i < n; ++i) f(i);
  end_code();
#endif
}

template <typename T>
typename T::value_type parallel_reduce(LO n, T f, std::string const& name = "") {
  typedef typename T::value_type VT;
  VT result;
  f.init(result);
#ifdef OMEGA_H_USE_KOKKOSCORE
  if (n > 0) Kokkos::parallel_reduce(name, policy(n), f, result);
#else
  begin_code(name);
  for (LO i = 0; i < n; ++i) f(i, result);
  end_code();
#endif
  return result;
}

template <typename T>
void parallel_scan(LO n, T f, std::string const& name = "") {
#ifdef OMEGA_H_USE_KOKKOSCORE
  if (n > 0) Kokkos::parallel_scan(policy(n), f, name);
#else
  typedef typename T::value_type VT;
  begin_code(name);
  VT update;
  f.init(update);
  for (LO i = 0; i < n; ++i) f(i, update, true);
  end_code();
#endif
}

}  // end namespace Omega_h

#endif
