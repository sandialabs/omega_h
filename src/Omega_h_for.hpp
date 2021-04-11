#ifndef OMEGA_H_FOR_HPP
#define OMEGA_H_FOR_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_profile.hpp>

#include <Omega_h_int_iterator.hpp>
#include <Omega_h_shared_alloc.hpp>

#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_kokkos.hpp>
#endif

#ifdef OMEGA_H_USE_CUDA

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <thrust/execution_policy.h>
#include <thrust/for_each.h>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#endif

namespace Omega_h {

#if defined(OMEGA_H_USE_CUDA)

namespace details {

template <class F, class ForwardIt>
__global__
void cuda_for_each(F f, ForwardIt first, ForwardIt last) {
  using difference_type = typename std::iterator_traits<ForwardIt>::difference_type;
  auto const i = static_cast<difference_type>(
          threadIdx.x + blockIdx.x * blockDim.x);
  ForwardIt const it = first + i;
  if (it < last) f(*it);
}

template <class T>
__host__ __device__ inline constexpr
T ceildiv(T a, T b) {
  return (a / b) + ((a % b) ? 1 : 0);
}

}

template <class ForwardIt, class UnaryFunction>
P3A_NEVER_INLINE
void for_each(
    ForwardIt first,
    ForwardIt last,
    UnaryFunction f)
{
  auto const n = last - first;
  if (n == 0) return;
  dim3 const cuda_block(32, 1, 1);
  dim3 const cuda_grid(details::ceildiv(unsigned(n), cuda_block.x), 1, 1);
  std::size_t const shared_memory_bytes = 0;
  cudaStream_t const cuda_stream = nullptr;
  details::cuda_for_each<<<
    cuda_grid,
    cuda_block,
    shared_memory_bytes,
    cuda_stream>>>(f, first, last);
}

#else

template <typename InputIterator, typename UnaryFunction>
void for_each(InputIterator first, InputIterator last, UnaryFunction&& f) {
  if (first >= last) return;
  Omega_h::entering_parallel = true;
  auto const f2 = std::move(f);
  Omega_h::entering_parallel = false;
#if defined(OMEGA_H_USE_CUDA)
  thrust::for_each(thrust::device, first, last, f2);
#elif defined(OMEGA_H_USE_OPENMP)
  LO const n = last - first;
#pragma omp parallel for
  for (LO i = 0; i < n; ++i) {
    f2(first[i]);
  }
#else
  for (; first != last; ++first) {
    f2(*first);
  }
#endif
}

#endif

template <typename UnaryFunction>
void parallel_for(LO n, UnaryFunction&& f) {
  auto const first = IntIterator(0);
  auto const last = IntIterator(n);
  ::Omega_h::for_each(first, last, f);
}

template <typename T>
void parallel_for(LO n, T const& f, char const* name = "") {
#if defined(OMEGA_H_USE_KOKKOS)
  if (n > 0) Kokkos::parallel_for(name, policy(n), f);
#else
  (void)name;
  auto f2 = f;
  parallel_for(n, std::move(f));
#endif
}

template <typename T>
void parallel_for(char const* name, LO n, T&& f) {
#if defined(OMEGA_H_USE_KOKKOS)
  if (n > 0) Kokkos::parallel_for(name, policy(n), f);
#else
  (void)name;
  parallel_for(n, std::move(f));
#endif
}

}  // end namespace Omega_h

#endif
