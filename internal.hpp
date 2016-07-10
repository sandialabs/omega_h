#ifndef INTERNAL_HPP
#define INTERNAL_HPP

#include "omega_h.hpp"

#define NORETURN(x) \
  do {              \
    assert(0);      \
    return x;       \
  } while (0)

#ifdef OSH_USE_CUDA
#define CONSTANT __constant__
#else
#define CONSTANT
#endif

#define CHECK(cond) OSH_CHECK(cond)
#define LAMBDA OSH_LAMBDA
#define INLINE OSH_INLINE
#define DEVICE OSH_DEVICE

#define EPSILON 1e-10
#define MANTISSA_BITS 52

#endif
