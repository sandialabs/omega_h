#ifndef OMEGA_H_REDUCE_HPP
#define OMEGA_H_REDUCE_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_stack.hpp>

#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_kokkos.hpp>
#elif defined(OMEGA_H_USE_CUDA)
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#endif

namespace Omega_h {

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

#if defined(OMEGA_H_USE_CUDA)

template <class Op>
Op native_op(Op const& op) const { return op; }
template<class T>
thrust::logical_and<T> native_op(Omega_h::logical_and<T> const&) const { return thrust::logical_and<T>(); }
template<class T>
thrust::plus<T> native_op(Omega_h::plus<T> const&) const { return thrust::plus<T>(); }
template <class T>
thrust::maximum<T> native_op()(Omega_h::maximum<T> const&) const { return thrust::maximum<T>(); }
template <class T>
thrust::minimum<T> native_op(Omega_h::minimum<T> const&) const { return thrust::minimum<T>(); }
template <class T>
thrust::identity<T> native_op(Omega_h::identity<T> const&) const { return thrust::identity<T>(); }

template <class Iterator, class Tranform, class Result, class Op>
Result transform_reduce(
    Iterator first,
    Iterator last,
    Tranform transform,
    Result init,
    Op op) {
  return thrust::transform_reduce(thrust::device, first, last, native_op(transform), init, native_op(op));
}

#elif defined(OMEGA_H_USE_OPENMP)

template <class Iterator, class Tranform, class Result>
Result transform_reduce(
    Iterator first,
    Iterator last,
    Tranform transform,
    Result init,
    Omega_h::logical_and<Result> op) {
  auto n = last - first;
#pragma omp parallel for reduction (&&:init)
  for (decltype(n) i = 0; i < n; ++i) {
    init = op(init, transform(first[i]));
  }
  return init;
}

template <class Iterator, class Tranform, class Result>
Result transform_reduce(
    Iterator first,
    Iterator last,
    Tranform transform,
    Result init,
    Omega_h::plus<Result> op) {
  auto n = last - first;
#pragma omp parallel for reduction (+:init)
  for (decltype(n) i = 0; i < n; ++i) {
    init = op(init, transform(first[i]));
  }
  return init;
}

template <class Iterator, class Tranform, class Result>
Result transform_reduce(
    Iterator first,
    Iterator last,
    Tranform transform,
    Result init,
    Omega_h::maximum<Result> op) {
  auto n = last - first;
#pragma omp parallel for reduction (max:init)
  for (decltype(n) i = 0; i < n; ++i) {
    init = op(init, transform(first[i]));
  }
  return init;
}

template <class Iterator, class Tranform, class Result>
Result transform_reduce(
    Iterator first,
    Iterator last,
    Tranform transform,
    Result init,
    Omega_h::minimum<Result> op) {
  auto n = last - first;
#pragma omp parallel for reduction (min:init)
  for (decltype(n) i = 0; i < n; ++i) {
    init = op(init, transform(first[i]));
  }
  return init;
}

#else

template <class Iterator, class Tranform, class Result, class Op>
Result transform_reduce(
    Iterator first,
    Iterator last,
    Tranform transform,
    Result init,
    Op op) {
  for (; first != last; ++first) {
    init = op(std::move(init), transform(*first));
  }
  return init;
}

#endif

}

#endif
