#ifndef OMEGA_H_KOKKOS_HPP
#define OMEGA_H_KOKKOS_HPP

#include <Omega_h_config.h>
#include <Omega_h_macros.h>
#include <Omega_h_defines.hpp>

OMEGA_H_SYSTEM_HEADER

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#if (__GNUC__ >= 7)
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wduplicated-branches"
#endif
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#if defined(KOKKOS_HAVE_CUDA) && !defined(OMEGA_H_USE_CUDA)
#error "Kokkos has CUDA, please reconfigure with Omega_h_USE_CUDA=ON"
#endif

namespace Omega_h {

#if defined(OMEGA_H_USE_CUDA)
  using ExecSpace = Kokkos::Cuda;
#elif defined(OMEGA_H_USE_HIP)
  using ExecSpace = Kokkos::HIP;
#elif defined(OMEGA_H_USE_SYCL)
  using ExecSpace = Kokkos::Experimental::SYCL;
#elif defined(OMEGA_H_USE_OpenMP)
  using ExecSpace = Kokkos::OpenMP;
#else
  using ExecSpace = Kokkos::DefaultExecutionSpace;
#endif

#if defined(OMEGA_H_MEM_SPACE_DEVICE)
  #if defined(OMEGA_H_USE_CUDA)
    using Space = Kokkos::CudaSpace;
  #elif defined(OMEGA_H_USE_HIP)
    using Space = Kokkos::HIPSpace;
  #elif defined(OMEGA_H_USE_SYCL)
    using Space = Kokkos::Experimental::SYCLDeviceUSMSpace;
  #endif
#elif defined(OMEGA_H_MEM_SPACE_SHARED)
  using Space = Kokkos::SharedSpace;
#elif defined(OMEGA_H_MEM_SPACE_HOSTPINNED)
  using Space = Kokkos::SharedHostPinnedSpace;
#else
  using Space = ExecSpace::memory_space;
#endif

using Device = Kokkos::Device<ExecSpace, Space>;
using StaticSched = Kokkos::Schedule<Kokkos::Static>;
using Policy = Kokkos::RangePolicy<ExecSpace, StaticSched, Omega_h::LO>;

template <class T> using View = Kokkos::View<T, Device>;

inline Policy policy(LO n) { return Policy(0, n); }
}  // namespace Omega_h

#endif
