#ifndef OMEGA_H_C_H
#define OMEGA_H_C_H

#include <assert.h>

#include <Omega_h_macros.h>
#include <Omega_h_mpi.h>

enum Omega_h_Type {
  OMEGA_H_I8 = 0,
  OMEGA_H_I32 = 2,
  OMEGA_H_I64 = 3,
  OMEGA_H_F64 = 5,
  OMEGA_H_REAL = OMEGA_H_F64,
};

enum { OMEGA_H_DIMS = 4 };

enum { OMEGA_H_VERT = 0, OMEGA_H_EDGE = 1, OMEGA_H_TRI = 2, OMEGA_H_TET = 3 };

enum Omega_h_Op { OMEGA_H_MIN, OMEGA_H_MAX, OMEGA_H_SUM };

enum Omega_h_Comparison { OMEGA_H_SAME, OMEGA_H_MORE, OMEGA_H_DIFF };

enum Omega_h_Transfer {
  OMEGA_H_INHERIT,
  OMEGA_H_LINEAR_INTERP,
  OMEGA_H_METRIC,
  OMEGA_H_CONSERVE,
  OMEGA_H_MOMENTUM_VELOCITY,
  OMEGA_H_POINTWISE,
};

enum Omega_h_Parting {
  OMEGA_H_ELEM_BASED,
  OMEGA_H_GHOSTED,
  OMEGA_H_VERT_BASED,
};

enum Omega_h_Source {
  OMEGA_H_HESSIAN,
  OMEGA_H_GIVEN,
  OMEGA_H_IMPLIED,
  OMEGA_H_PROXIMITY,
  OMEGA_H_CURVATURE,
};

enum Omega_h_Scales {
  OMEGA_H_ABSOLUTE,
  OMEGA_H_SCALES,
};

#ifdef __cplusplus
extern "C" {
#endif

void Omega_h_protect();
void Omega_h_signal_handler(int s);
void Omega_h_fail(char const* format, ...)
    __attribute__((noreturn, format(printf, 1, 2)));

#ifdef __cplusplus
}  // end of extern "C" block
#endif

#ifdef __CUDA_ARCH__
#define OMEGA_H_CHECK(cond) assert(cond)
#else
#define OMEGA_H_CHECK(cond)                                                    \
  ((cond) ? ((void)0)                                                          \
          : Omega_h_fail(                                                      \
                "assertion %s failed at %s +%d\n", #cond, __FILE__, __LINE__))
#endif

#ifdef __clang__
#define OMEGA_H_NORETURN(x) assert(false)
#else
#define OMEGA_H_NORETURN(x)                                                    \
  do {                                                                         \
    assert(false);                                                             \
    return x;                                                                  \
  } while (false)
#endif

#endif
