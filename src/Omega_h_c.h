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
};

enum { OMEGA_H_DIMS = 4 };

enum { OMEGA_H_VERT = 0, OMEGA_H_EDGE = 1, OMEGA_H_TRI = 2, OMEGA_H_TET = 3 };

enum Omega_h_Op { OMEGA_H_MIN, OMEGA_H_MAX, OMEGA_H_SUM, OMEGA_H_BIT_OR };

enum Omega_h_Xfer {
  OMEGA_H_DONT_TRANSFER      = 0,
  OMEGA_H_INHERIT            = 1,
  OMEGA_H_LINEAR_INTERP      = 2,
  OMEGA_H_POINTWISE          = 3,
  OMEGA_H_CONSERVE           = 4,
  OMEGA_H_GLOBAL             = 5,
  OMEGA_H_LENGTH             = 6,
  OMEGA_H_QUALITY            = 7,
  OMEGA_H_METRIC             = 8,
/*OMEGA_H_CONSERVE_R3D       = 9,*/
  OMEGA_H_MOMENTUM_VELOCITY  =10,
  OMEGA_H_SIZE               =11
};

enum Omega_h_Parting {
  OMEGA_H_ELEM_BASED,
  OMEGA_H_GHOSTED,
  OMEGA_H_VERT_BASED,
};

enum Omega_h_Comparison { OMEGA_H_SAME, OMEGA_H_MORE, OMEGA_H_DIFF };

enum Omega_h_Outflags {
  OMEGA_H_DONT_OUTPUT = 0x0,
  OMEGA_H_DO_SAVE = 0x1,
  OMEGA_H_DO_VIZ = 0x2,
  OMEGA_H_DO_OUTPUT = OMEGA_H_DO_SAVE | OMEGA_H_DO_VIZ
};

#ifdef __cplusplus
extern "C" {
#endif

void Omega_h_fail(char const* format, ...)
    __attribute__((noreturn, format(printf, 1, 2)));

#ifdef __cplusplus
}  // end of extern "C" block
#endif

#ifdef __CUDA_ARCH__
#define OMEGA_H_CHECK(cond) assert(cond)
#else
#define OMEGA_H_CHECK(cond)                                                    \
  ((cond) ? ((void)0) : Omega_h_fail("assertion %s failed at %s +%d\n", #cond, \
                            __FILE__, __LINE__))
#endif

#endif
