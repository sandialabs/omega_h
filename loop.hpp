template <typename T>
void parallel_for(Int n, T const& f) {
#ifdef USE_KOKKOS
  Kokkos::parallel_for(n, f);
#else
  for (Int i = 0; i < n; ++i)
    f(i);
#endif
}

template <typename T>
typename T::value_type parallel_reduce(Int n, T f)
{
  typename T::value_type result;
  f.init(result);
#ifdef USE_KOKKOS
  Kokkos::parallel_reduce(n, f, result);
#else
  for (Int i = 0; i < n; ++i)
    f(i, result);
#endif
  return result;
}

template <typename T>
void parallel_scan(Int n, T f)
{
#ifdef USE_KOKKOS
  Kokkos::parallel_scan(n, f);
#else
  typename T::value_type update;
  f.init(update);
  for (Int i = 0; i < n; ++i)
    f(i, update, true);
#endif
}
