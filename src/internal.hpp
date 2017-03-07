#ifndef INTERNAL_HPP
#define INTERNAL_HPP

#include "Omega_h.hpp"

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

/* The environment on Sandia's `curie` machine won't
 * give us this through <cstdint> even though we're
 * in C++11 mode.
 */
#ifndef INT8_MAX
#define INT8_MAX 0x7f
#endif

#endif
