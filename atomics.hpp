template <class T>
OSH_INLINE void atomic_increment(volatile T* const dest)
{
#ifdef OSH_USE_KOKKOS
  return Kokkos::atomic_increment(dest);
#else
  ++(*dest);
#endif
}

template <class T>
OSH_INLINE void atomic_add(volatile T* const dest, const T val)
{
#ifdef OSH_USE_KOKKOS
  return Kokkos::atomic_add(dest, val);
#else
  *dest += val;
#endif
}

template <class T>
OSH_INLINE T atomic_fetch_add(volatile T* const dest, const T val)
{
#ifdef OSH_USE_KOKKOS
  return Kokkos::atomic_fetch_add(dest, val);
#else
  T tmp = *dest;
  *dest += val;
  return tmp;
#endif
}
