#ifndef OMEGA_H_REDUCE_HPP
#define OMEGA_H_REDUCE_HPP

#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_kokkos.hpp>
#endif

#include <Omega_h_scalar.hpp>
#include <Omega_h_shared_alloc.hpp>
#if defined(OMEGA_H_USE_CUDA)
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#endif

namespace Omega_h {

template <typename T>
typename T::value_type parallel_reduce(LO n, T f, char const* name = "") {
  using VT = typename T::value_type;
#ifdef OMEGA_H_USE_KOKKOS
  VT result;
  f.init(result);
  if (n > 0) Kokkos::parallel_reduce(name, policy(n), f, result);
#else
  (void)name;
  VT result;
  f.init(result);
  for (LO i = 0; i < n; ++i) f(i, result);
#endif
  return result;
}

#if defined(OMEGA_H_USE_CUDA)

template <class Op>
Op native_op(Op const& op) {
  return op;
}
template <class T>
thrust::logical_and<T> native_op(Omega_h::logical_and<T> const&) {
  return thrust::logical_and<T>();
}
template <class T>
thrust::plus<T> native_op(Omega_h::plus<T> const&) {
  return thrust::plus<T>();
}
template <class T>
thrust::maximum<T> native_op(Omega_h::maximum<T> const&) {
  return thrust::maximum<T>();
}
template <class T>
thrust::minimum<T> native_op(Omega_h::minimum<T> const&) {
  return thrust::minimum<T>();
}
template <class T>
thrust::identity<T> native_op(Omega_h::identity<T> const&) {
  return thrust::identity<T>();
}

template <class Iterator, class Tranform, class Result, class Op>
Result transform_reduce(
    Iterator first, Iterator last, Result init, Op op, Tranform transform) {
  Omega_h::entering_parallel = true;
  auto const transform_parallel = std::move(transform);
  Omega_h::entering_parallel = false;
  return thrust::transform_reduce(thrust::device, first, last,
      native_op(transform_parallel), init, native_op(op));
}

#elif defined(OMEGA_H_USE_OPENMP)

template <class Iterator, class Tranform, class Result, class Op>
Result transform_reduce(
    Iterator first, Iterator last, Result init, Op, Tranform&& transform) {
#pragma omp declare reduction(osh:Result                                       \
                              : omp_out = Op()(omp_out, omp_in))               \
    initializer(omp_priv = omp_orig)
  LO const n = last - first;
  Omega_h::entering_parallel = true;
  auto const transform_local = std::move(transform);
  Omega_h::entering_parallel = false;
#pragma omp parallel for reduction(osh : init)
  for (LO i = 0; i < n; ++i) {
    init = Op()(init, transform_local(first[i]));
  }
  return init;
}

template <class Iterator, class Tranform, class Result>
Result transform_reduce(Iterator first, Iterator last, Result init,
    Omega_h::logical_and<Result> op, Tranform&& transform) {
  LO const n = last - first;
  Omega_h::entering_parallel = true;
  auto const transform_local = std::move(transform);
  Omega_h::entering_parallel = false;
#pragma omp parallel for reduction(&& : init)
  for (LO i = 0; i < n; ++i) {
    init = op(init, transform_local(first[i]));
  }
  return init;
}

template <class Iterator, class Tranform, class Result>
Result transform_reduce(Iterator first, Iterator last, Result init,
    Omega_h::plus<Result> op, Tranform&& transform) {
  LO const n = LO(last - first);
  Omega_h::entering_parallel = true;
  auto const transform_local = std::move(transform);
  Omega_h::entering_parallel = false;
#pragma omp parallel for reduction(+ : init)
  for (LO i = 0; i < n; ++i) {
    init = op(init, transform_local(first[i]));
  }
  return init;
}

template <class Iterator, class Tranform, class Result>
Result transform_reduce(Iterator first, Iterator last, Result init,
    Omega_h::maximum<Result> op, Tranform&& transform) {
  LO const n = last - first;
  Omega_h::entering_parallel = true;
  auto const transform_local = std::move(transform);
  Omega_h::entering_parallel = false;
#pragma omp parallel for reduction(max : init)
  for (LO i = 0; i < n; ++i) {
    init = op(init, transform_local(first[i]));
  }
  return init;
}

template <class Iterator, class Tranform, class Result>
Result transform_reduce(Iterator first, Iterator last, Result init,
    Omega_h::minimum<Result> op, Tranform&& transform) {
  LO const n = last - first;
  Omega_h::entering_parallel = true;
  auto const transform_local = std::move(transform);
  Omega_h::entering_parallel = false;
#pragma omp parallel for reduction(min : init)
  for (LO i = 0; i < n; ++i) {
    init = op(init, transform_local(first[i]));
  }
  return init;
}

#else

template <class Iterator, class Tranform, class Result, class Op>
Result transform_reduce(
    Iterator first, Iterator last, Result init, Op op, Tranform&& transform) {
  Omega_h::entering_parallel = true;
  auto const transform_local = std::move(transform);
  Omega_h::entering_parallel = false;
  for (; first != last; ++first) {
    init = op(std::move(init), transform_local(*first));
  }
  return init;
}

#endif

}  // namespace Omega_h

#endif
