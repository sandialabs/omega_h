#ifndef OMEGA_H_FOR_HPP
#define OMEGA_H_FOR_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_stack.hpp>

#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_kokkos.hpp>
#elif defined(OMEGA_H_USE_OPENMP) || defined(OMEGA_H_USE_CUDA)
#include <Omega_h_shared_alloc.hpp>
#endif

namespace Omega_h {

#if defined(OMEGA_H_USE_CUDA) && (!defined(OMEGA_H_USE_KOKKOSCORE))
template <typename T>
__global__ static void launch_cuda(const T f, const LO n) {
  LO i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  f(i);
}
extern int block_size_cuda;
#endif

template <typename T>
void parallel_for(LO n, T const& f, char const* name = "") {
#if defined(OMEGA_H_USE_KOKKOSCORE)
  if (n > 0) Kokkos::parallel_for(name, policy(n), f);
#elif defined(OMEGA_H_USE_CUDA)
  begin_code(name);
  entering_parallel = true;
  T f2 = f;
  entering_parallel = false;
  LO nblocks = (n + block_size_cuda - 1) / block_size_cuda;
  launch_cuda<T><<<nblocks, block_size>>>(f2, n);
  end_code();
#elif defined(OMEGA_H_USE_OPENMP)
  begin_code(name);
  entering_parallel = true;
  const T f2 = f;
  entering_parallel = false;
#pragma omp parallel for
  for (LO i = 0; i < n; ++i) {
    f2(i);
  }
  end_code();
#else
  begin_code(name);
  const T& f2 = f;
  for (LO i = 0; i < n; ++i) f2(i);
  end_code();
#endif
}

template <typename T>
void parallel_for(char const* name, LO n, T&& f) {
#if defined(OMEGA_H_USE_KOKKOSCORE)
  if (n > 0) Kokkos::parallel_for(name, policy(n), f);
#elif defined(OMEGA_H_USE_CUDA)
  begin_code(name);
  entering_parallel = true;
  T f2 = std::move(f);
  entering_parallel = false;
  LO nblocks = (n + block_size_cuda - 1) / block_size_cuda;
  launch_cuda<T><<<nblocks, block_size>>>(f2, n);
  end_code();
#elif defined(OMEGA_H_USE_OPENMP)
  begin_code(name);
  entering_parallel = true;
  const T f2 = std::move(f);
  entering_parallel = false;
#pragma omp parallel for
  for (LO i = 0; i < n; ++i) {
    f2(i);
  }
  end_code();
#else
  begin_code(name);
  const T& f2 = f;
  for (LO i = 0; i < n; ++i) f2(i);
  end_code();
#endif
}

// this class exists simply to allow users to omit a closing parenthesis:
// OMEGA_H_FOR("name", i, n) {
// }; // <- no closing paren needed!
// as opposed to
// parallel_for("name", n, OMEGA_H_LAMBDA(LO i) {
// }); // <- closing paren needed!
struct ParallelFor {
  char const* name;
  LO n;
  ParallelFor(char const* name_in, LO n_in):name(name_in),n(n_in) {}
  template <typename T>
  void operator<<(T&& f) {
    parallel_for(name, n, std::move(f));
  }
};

#define OMEGA_H_FOR(name, i, n) ParallelFor(name, n) << OMEGA_H_LAMBDA(Omega_h::LO i)

}  // end namespace Omega_h

#endif
