#ifndef OMEGA_H_FOR_HPP
#define OMEGA_H_FOR_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_stack.hpp>

#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_kokkos.hpp>
#else
#include <Omega_h_shared_alloc.hpp>
#endif

namespace Omega_h {

#if defined(OMEGA_H_USE_CUDA)
template <typename T>
__global__ static void launch_cuda(const T f, const LO n) {
  LO i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  f(i);
}
extern int block_size_cuda;
#endif

template <typename T>
void parallel_for(LO n, T&& f) {
#if defined(OMEGA_H_USE_CUDA)
  Omega_h::entering_parallel = true;
  auto const f2 = std::move(f);
  Omega_h::entering_parallel = false;
  LO nblocks = (n + block_size_cuda - 1) / block_size_cuda;
  launch_cuda<T><<<nblocks, block_size> > >(f2, n);
#elif defined(OMEGA_H_USE_OPENMP)
  Omega_h::entering_parallel = true;
  auto const f2 = std::move(f);
  Omega_h::entering_parallel = false;
#pragma omp parallel for
  for (LO i = 0; i < n; ++i) {
    f2(i);
  }
#else
  Omega_h::entering_parallel = true;
  auto const f2 = std::move(f);
  Omega_h::entering_parallel = false;
  for (LO i = 0; i < n; ++i) {
    f2(i);
  }
#endif
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
