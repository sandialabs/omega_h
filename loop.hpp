template <typename T>
void parallel_for(UInt n, T const& f) {
#ifdef USE_KOKKOS
  Kokkos::parallel_for(n, f);
#else
  for (UInt i = 0; i < n; ++i)
    f(i);
#endif
}

template <typename T>
typename T::value_type parallel_reduce(UInt n, T f)
{
  typename T::value_type result;
  f.init(result);
#ifdef USE_KOKKOS
  Kokkos::parallel_reduce(n, f, result);
#else
  for (UInt i = 0; i < n; ++i)
    f(i, result);
#endif
  return result;
}

template <typename T>
void parallel_scan(UInt n, T f)
{
#ifdef USE_KOKKOS
  Kokkos::parallel_scan(n, f);
#else
  typename T::value_type update;
  f.init(update);
  for (UInt i = 0; i < n; ++i)
    f(i, update, true);
#endif
}
