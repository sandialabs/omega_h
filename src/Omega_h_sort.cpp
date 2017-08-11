#include "Omega_h_sort.hpp"

#include <algorithm>

#if defined(OMEGA_H_USE_CUDA)
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#elif defined(OMEGA_H_USE_OPENMP)
#include <omp.h>
#include "intel_sort/parallel_stable_sort.hpp"
#include "intel_sort/pss_common.hpp"
#endif

#include "Omega_h_control.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_scalar.hpp"

namespace Omega_h {

template <typename T, typename Comp>
static void parallel_sort(T* b, T* e, Comp c) {
  auto t0 = now();
#if defined(OMEGA_H_USE_CUDA)
  auto bptr = thrust::device_ptr<T>(b);
  auto eptr = thrust::device_ptr<T>(e);
  thrust::stable_sort(bptr, eptr, c);
#elif defined(OMEGA_H_USE_OPENMP)
  pss::parallel_stable_sort(b, e, c);
#else
  std::stable_sort(b, e, c);
#endif
  auto t1 = now();
  add_to_global_timer("sorting", t1 - t0);
}

template <typename T, Int N>
struct CompareKeySets {
  T const* keys_;
  CompareKeySets(T const* keys) : keys_(keys) {}
  OMEGA_H_INLINE bool operator()(const LO& a, const LO& b) const {
    for (Int i = 0; i < N; ++i) {
      T x = keys_[a * N + i];
      T y = keys_[b * N + i];
      if (x != y) return x < y;
    }
    return false;
  }
};

template <Int N, typename T>
static LOs sort_by_keys_tmpl(Read<T> keys) {
  auto n = divide_no_remainder(keys.size(), N);
  Write<LO> perm(n, 0, 1);
  LO* begin = perm.data();
  LO* end = perm.data() + n;
  T const* keyptr = keys.data();
  parallel_sort<LO, CompareKeySets<T, N>>(
      begin, end, CompareKeySets<T, N>(keyptr));
  return perm;
}

template <typename T>
LOs sort_by_keys(Read<T> keys, Int width) {
  switch (width) {
    case 1:
      return sort_by_keys_tmpl<1>(keys);
    case 2:
      return sort_by_keys_tmpl<2>(keys);
    case 3:
      return sort_by_keys_tmpl<3>(keys);
  }
  OMEGA_H_NORETURN(LOs());
}

#define INST(T) template LOs sort_by_keys(Read<T> keys, Int width);
INST(LO)
INST(GO)
#undef INST

}  // end namespace Omega_h
