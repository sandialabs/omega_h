#ifndef OMEGA_H_MACROS_H
#define OMEGA_H_MACROS_H

#include <Omega_h_config.h>

#define OMEGA_H_STRINGIFY(s) #s
/* apparently you need two macros to make a string */
#define OMEGA_H_TOSTRING(s) OMEGA_H_STRINGIFY(s)

#define OMEGA_H_PRAGMA(x) _Pragma(#x)

#if defined(__clang__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(clang system_header)
#elif defined(__GNUC__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(GCC system_header)
#endif

#define OMEGA_H_PI 3.14159265358979323
//                  ^    ^    ^    ^
//                  0    5   10   15
#define OMEGA_H_EPSILON 1e-10
#define OMEGA_H_MANTISSA_BITS 52

#if defined(__clang__)
#define OMEGA_H_FALLTHROUGH [[clang::fallthrough]]
/* Test for GCC >= 7.0.0 */
#elif defined(__GNUC__) && (__GNUC__ > 7 || __GNUC__ == 7)
#define OMEGA_H_FALLTHROUGH [[gnu::fallthrough]]
#else
#define OMEGA_H_FALLTHROUGH ((void)0)
#endif

#define OMEGA_H_NODISCARD __attribute__((warn_unused_result))

#ifdef OMEGA_H_USE_CUDA
#define OMEGA_H_INLINE __host__ __device__ inline
#define OMEGA_H_DEVICE __device__ inline
#define OMEGA_H_LAMBDA [=] __device__
#define OMEGA_H_CONSTANT_DATA __constant__
#else
#define OMEGA_H_INLINE inline
#define OMEGA_H_DEVICE inline
#define OMEGA_H_LAMBDA [=]
#define OMEGA_H_CONSTANT_DATA
#endif

#endif
