#ifndef LOOP_HPP
#define LOOP_HPP

namespace osh {

template <typename T>
void parallel_for(Int n, T const& f) {
#ifdef OSH_USE_KOKKOS
  if (n > 0) Kokkos::parallel_for(static_cast<std::size_t>(n), f);
#else
  for (Int i = 0; i < n; ++i) f(i);
#endif
}

template <typename T>
typename T::value_type parallel_reduce(Int n, T f) {
  typename T::value_type result;
  f.init(result);
#ifdef OSH_USE_KOKKOS
  if (n > 0) Kokkos::parallel_reduce(static_cast<std::size_t>(n), f, result);
#else
  for (Int i = 0; i < n; ++i) f(i, result);
#endif
  return result;
}

template <typename T>
void parallel_scan(Int n, T f) {
#ifdef OSH_USE_KOKKOS
  if (n > 0) Kokkos::parallel_scan(static_cast<std::size_t>(n), f);
#else
  typename T::value_type update;
  f.init(update);
  for (Int i = 0; i < n; ++i) f(i, update, true);
#endif
}

} //end namespace osh

#endif
