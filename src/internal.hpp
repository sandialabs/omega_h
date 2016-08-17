#ifndef INTERNAL_HPP
#define INTERNAL_HPP

#include "omega_h.hpp"

#ifdef __clang__
#define NORETURN(x) assert(false)
#else
#define NORETURN(x)                                                            \
  do {                                                                         \
    assert(false);                                                             \
    return x;                                                                  \
  } while (false)
#endif

#ifdef OMEGA_H_USE_CUDA
#define CONSTANT __constant__
#else
#define CONSTANT
#endif

#define CHECK(cond) OMEGA_H_CHECK(cond)
#define LAMBDA OMEGA_H_LAMBDA
#define INLINE OMEGA_H_INLINE
#define DEVICE OMEGA_H_DEVICE

#define EPSILON 1e-10
#define MANTISSA_BITS 52

#endif
