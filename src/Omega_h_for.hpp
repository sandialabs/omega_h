#ifndef OMEGA_H_FOR_HPP
#define OMEGA_H_FOR_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_stack.hpp>

#include <Omega_h_shared_alloc.hpp>
#include <Omega_h_int_iterator.hpp>
#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_kokkos.hpp>
#else
#ifdef OMEGA_H_USE_CUDA
#include <thrust/for_each.h>
#include <thrust/execution_policy.h>
#endif
#endif

namespace Omega_h {

template <typename InputIterator, typename UnaryFunction>
void for_each(InputIterator first, InputIterator last, UnaryFunction&& f) {
  if (first >= last) return;
#if defined(OMEGA_H_USE_CUDA)
  Omega_h::entering_parallel = true;
  auto const f2 = std::move(f);
  Omega_h::entering_parallel = false;
  auto const first = IntIterator(0);
  auto const last = IntIterator(n);
  thrust::for_each(thrust::device, first, last, f2);
#elif defined(OMEGA_H_USE_OPENMP)
  LO const n = last - first;
  Omega_h::entering_parallel = true;
  auto const f2 = std::move(f);
  Omega_h::entering_parallel = false;
#pragma omp parallel for
  for (LO i = 0; i < n; ++i) {
    f2(first[i]);
  }
#else
  Omega_h::entering_parallel = true;
  auto const f2 = std::move(f);
  Omega_h::entering_parallel = false;
  for (; first != last; ++first) {
    f2(*first);
  }
#endif
}

template <typename UnaryFunction>
void parallel_for(LO n, UnaryFunction&& f) {
  auto const first = IntIterator(0);
  auto const last = IntIterator(n);
  for_each(first, last, f);
}

template <typename T>
void parallel_for(LO n, T const& f, char const* name = "") {
#if defined(OMEGA_H_USE_KOKKOSCORE)
  if (n > 0) Kokkos::parallel_for(name, policy(n), f);
#else
  (void)name;
  auto f2 = f;
  parallel_for(n, std::move(f));
#endif
}

template <typename T>
void parallel_for(char const* name, LO n, T&& f) {
#if defined(OMEGA_H_USE_KOKKOSCORE)
  if (n > 0) Kokkos::parallel_for(name, policy(n), f);
#else
  (void)name;
  parallel_for(n, std::move(f));
#endif
}

}  // end namespace Omega_h

#endif
