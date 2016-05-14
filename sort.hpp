#if defined(USE_CUDA)
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#elif defined(USE_OPENMP)
#include "bundled/parallel_stable_sort.hpp"
#else
#include <algorithm>
#endif

template <typename T, typename Comp>
void parallel_sort(T* b, T* e, Comp c) {
#if defined(USE_CUDA)
  thrust::stable_sort(
      thrust::device_ptr<T>(b),
      thrust::device_ptr<T>(e), c);
#elif defined(USE_OPENMP)
  pss::parallel_stable_sort(b, e, c);
#else
  std::stable_sort(b, e, c);
#endif
}

template <typename T>
struct Less {
  bool operator()(T const& a, T const& b) {
    return a < b;
  }
};

template <typename T, typename Comp>
void parallel_sort(T* b, T* e) {
  parallel_sort(b, e, Less<T>());
}
