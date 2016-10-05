#ifndef OMEGA_H_MACROS_H
#define OMEGA_H_MACROS_H

#include <Omega_h_config.h>

/* this block of preprocessor code creates a string
   literal describing the version and compile options
   used to configure this header.
   this is used to ensure that during user application
   compiles, the header included matches the library linked to.
   it is not foolproof, but better than nothing */
#define OMEGA_H_TOSTR2(s) #s
#define OMEGA_H_TOSTR(s) OMEGA_H_TOSTR2(s)
#ifdef OMEGA_H_USE_MPI
#define OMEGA_H_MPI_STR "1"
#else
#define OMEGA_H_MPI_STR "0"
#endif
#ifdef OMEGA_H_USE_KOKKOS
#define OMEGA_H_KOKKOS_STR "1"
#else
#define OMEGA_H_KOKKOS_STR "0"
#endif
#ifdef OMEGA_H_USE_OPENMP
#define OMEGA_H_OPENMP_STR "1"
#else
#define OMEGA_H_OPENMP_STR "0"
#endif
#ifdef OMEGA_H_USE_CUDA
#define OMEGA_H_CUDA_STR "1"
#else
#define OMEGA_H_CUDA_STR "0"
#endif
#ifdef OMEGA_H_USE_ZLIB
#define OMEGA_H_ZLIB_STR "1"
#else
#define OMEGA_H_ZLIB_STR "0"
#endif
#ifdef OMEGA_H_CHECK_BOUNDS
#define OMEGA_H_BOUNDS_STR "1"
#else
#define OMEGA_H_BOUNDS_STR "0"
#endif
#define OMEGA_H_VERSION                                 \
  OMEGA_H_TOSTR(OMEGA_H_VERSION_MAJOR) "."              \
  OMEGA_H_TOSTR(OMEGA_H_VERSION_MINOR) "."              \
  OMEGA_H_TOSTR(OMEGA_H_VERSION_PATCH) "+"              \
  OMEGA_H_MPI_STR OMEGA_H_KOKKOS_STR OMEGA_H_OPENMP_STR \
  OMEGA_H_CUDA_STR OMEGA_H_ZLIB_STR OMEGA_H_BOUNDS_STR

#define OMEGA_H_PRAGMA(x) _Pragma(#x)

#if defined(__clang__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(clang system_header)
#elif defined(__GNUC__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(GCC system_header)
#endif

#endif
