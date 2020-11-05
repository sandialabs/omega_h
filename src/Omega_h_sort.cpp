#include <Omega_h_int_iterator.hpp>
#include <Omega_h_scan.hpp>
#include <Omega_h_sort.hpp>

#include <algorithm>
#include <vector>

#if defined(OMEGA_H_USE_CUDA)

#if defined(__clang__)
template <class... Args>
inline __host__ __device__ void va_printf(const char*, Args...) {
  printf("\n");
}
#endif
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#elif defined(OMEGA_H_USE_OPENMP)

#include <omp.h>
#include <pss/parallel_stable_sort.hpp>
#include <pss/pss_common.hpp>

#endif

#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_scalar.hpp"
#include "Omega_h_timer.hpp"

namespace Omega_h {

template <typename T, typename Comp>
static void parallel_sort(T* b, T* e, Comp c) {
  begin_code("parallel_sort");
#if defined(OMEGA_H_USE_CUDA)
  auto bptr = thrust::device_ptr<T>(b);
  auto eptr = thrust::device_ptr<T>(e);
  thrust::stable_sort(bptr, eptr, c);
#elif defined(OMEGA_H_USE_OPENMP)
  pss::parallel_stable_sort(b, e, c);
#else
  std::stable_sort(b, e, c);
#endif
  end_code();
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
  begin_code("sort_by_keys");
  auto n = divide_no_remainder(keys.size(), N);
  Write<LO> perm(n, 0, 1);
  LO* begin = perm.data();
  LO* end = perm.data() + n;
  T const* keyptr = keys.data();
  parallel_sort<LO, CompareKeySets<T, N>>(
      begin, end, CompareKeySets<T, N>(keyptr));
  end_code();
  return perm;
}

template <typename T>
LOs sort_by_keys(Read<T> keys, Int width) {
  if (width == 1) return sort_by_keys_tmpl<1>(keys);
  if (width == 2) return sort_by_keys_tmpl<2>(keys);
  if (width == 3) return sort_by_keys_tmpl<3>(keys);
  if (width == 4) return sort_by_keys_tmpl<4>(keys);
  OMEGA_H_NORETURN(LOs());
}

#define INST(T) template LOs sort_by_keys(Read<T> keys, Int width);
INST(LO)
INST(GO)
#undef INST

template <typename T>
T next_smallest_value(Read<T> const a, T const value) {
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const init = ArithTraits<T>::max();
  auto const op = minimum<T>();
  auto transform = OMEGA_H_LAMBDA(LO i)->T {
    return (a[i] > value) ? a[i] : init;
  };
  return transform_reduce(first, last, init, op, std::move(transform));
}

template <typename T>
LO number_same_values(
    Read<T> const a, T const value, Write<LO> const tmp_perm) {
  tmp_perm.set(0, 0);
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const result = tmp_perm.begin() + 1;
  auto const op = plus<LO>();
  auto transform = OMEGA_H_LAMBDA(LO i)->LO {
    return a[i] == value ? LO(1) : LO(0);
  };
  transform_inclusive_scan(first, last, result, op, std::move(transform));
  return read(tmp_perm).last();
}

template <typename T>
void combine_perms(Read<T> const a, T const value,
    Write<LO> const tmp_perm, Write<LO> const perm, LO ndone) {
  tmp_perm.set(0, 0);
  auto functor = OMEGA_H_LAMBDA(LO i) {
    if (a[i] == value) perm[i] = tmp_perm[i] + ndone;
  };
  parallel_for(a.size(), std::move(functor));
}

template <typename T>
void sort_small_range(Read<T> a, LOs* p_perm, LOs* p_fan, Read<T>* p_uniq) {
  LO ndone = 0;
  T value = get_min(a);
  Write<LO> tmp_perm(a.size() + 1);
  Write<LO> perm(a.size(), -42);
  std::vector<LO> fan_vec;
  std::vector<T> uniq_vec;
  fan_vec.push_back(0);
  if (a.size()) {
    while (true) {
      uniq_vec.push_back(value);
      auto const ndid = number_same_values(a, value, tmp_perm);
      combine_perms(a, value, tmp_perm, perm, ndone);
      ndone += ndid;
      fan_vec.push_back(ndone);
      if (ndone == a.size()) break;
      value = next_smallest_value(a, value);
    }
  }
  HostWrite<LO> h_fan(LO(fan_vec.size()));
  for (LO i = 0; i < h_fan.size(); ++i) {
    h_fan[i] = fan_vec[std::size_t(i)];
  }
  HostWrite<T> h_uniq(LO(uniq_vec.size()));
  for (LO i = 0; i < h_uniq.size(); ++i) {
    h_uniq[i] = uniq_vec[std::size_t(i)];
  }
  *p_perm = perm;
  *p_fan = h_fan.write();
  *p_uniq = h_uniq.write();
}

template void sort_small_range(
    Read<I32> items2values, LOs* p_perm, LOs* p_fan, Read<I32>* p_uniq);

}  // end namespace Omega_h
