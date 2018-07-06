#include "Omega_h_sort.hpp"

#include <algorithm>
#include <vector>

#if defined(OMEGA_H_USE_CUDA)
#if defined(__clang__)
template <class... Args>
inline __host__ __device__ void va_printf(const char*, Args...) {
  printf("\n");
}
#endif
#include <Omega_h_thrust.hpp>
#elif defined(OMEGA_H_USE_OPENMP)
#include <omp.h>
#include <pss/parallel_stable_sort.hpp>
#include <pss/pss_common.hpp>
#endif

#include "Omega_h_array_ops.hpp"
#include "Omega_h_loop.hpp"
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
struct SortSmallRangeValue {
  T next_smallest;
  LO index;
};

template <typename T>
struct SortSmallRange {
  using value_type = SortSmallRangeValue<T>;
  using input_type = T;
  Read<T> items2values_;
  Write<LO> perm_;
  LO nitems_;
  LO target_value_;
  Write<LO> ndone_after_;
  Write<T> next_smallest_;
  LO ndone_before_;
  OMEGA_H_INLINE void init(value_type& update) const {
    update.next_smallest = ArithTraits<T>::max();
    update.index = 0;
  }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update.next_smallest = min2(update.next_smallest, input.next_smallest);
    update.index += input.index;
  }
  OMEGA_H_DEVICE void operator()(
      LO i, value_type& update, bool final_pass) const {
    auto value = items2values_[i];
    if (value == target_value_) {
      ++update.index;
    } else if (value > target_value_) {
      update.next_smallest = min2(update.next_smallest, value);
    }
    if (final_pass) {
      if (i == nitems_ - 1) {
        next_smallest_[0] = update.next_smallest;
        ndone_after_[0] = update.index + ndone_before_;
      }
      if (value == target_value_) {
        perm_[i] = update.index + ndone_before_ - 1;
      }
    }
  }
  void run(Read<T> items2values, LOs* p_perm, LOs* p_fan, Read<T>* p_uniq) {
    nitems_ = items2values.size();
    if (nitems_ == 0) {
      *p_perm = LOs({});
      *p_fan = LOs({0});
      *p_uniq = LOs({});
      return;
    }
    items2values_ = items2values;
    ndone_after_ = Write<LO>(1, 0);
    next_smallest_ = Write<LO>(1);
    perm_ = Write<LO>(nitems_);
    target_value_ = get_min(items2values);
    ndone_before_ = 0;
    std::vector<T> uniq_vector;
    std::vector<T> fan_vector;
    while (ndone_before_ < nitems_) {
      uniq_vector.push_back(target_value_);
      fan_vector.push_back(ndone_before_);
      ::Omega_h::parallel_scan(nitems_, *this, "sort_small_range");
      OMEGA_H_CHECK(ndone_after_.get(0) > ndone_before_);
      target_value_ = next_smallest_.get(0);
      ndone_before_ = ndone_after_.get(0);
    }
    fan_vector.push_back(ndone_before_);
    auto nuniq = LO(uniq_vector.size());
    auto nfan = LO(fan_vector.size());
    OMEGA_H_CHECK(nfan == nuniq + 1);
    auto uniq_h = HostWrite<LO>(nuniq);
    auto fan_h = HostWrite<LO>(nfan);
    for (LO i = 0; i < nuniq; ++i) {
      uniq_h[i] = uniq_vector[std::size_t(i)];
      fan_h[i] = fan_vector[std::size_t(i)];
    }
    fan_h[nuniq] = fan_vector[std::size_t(nuniq)];
    *p_perm = perm_;
    *p_uniq = uniq_h.write();
    *p_fan = fan_h.write();
  }
};

template <typename T>
void sort_small_range(
    Read<T> items2values, LOs* p_perm, LOs* p_fan, Read<T>* p_uniq) {
  SortSmallRange<T> f;
  f.run(items2values, p_perm, p_fan, p_uniq);
}

template void sort_small_range(
    Read<I32> items2values, LOs* p_perm, LOs* p_fan, Read<I32>* p_uniq);

}  // end namespace Omega_h
