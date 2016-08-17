#ifndef OMEGA_H_C_H
#define OMEGA_H_C_H

#include <assert.h>

#include "Omega_h_config.h"

#define OMEGA_H_PRAGMA(x) _Pragma(#x)

#if defined(__clang__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(clang system_header)
#elif defined(__GNUC__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(GCC system_header)
#endif

#include "Omega_h_mpi.h"

enum Omega_h_Type {
  OMEGA_H_I8 = 0,
  OMEGA_H_I32 = 2,
  OMEGA_H_I64 = 3,
  OMEGA_H_F64 = 5,
};

enum { OMEGA_H_DIMS = 4 };

enum { OMEGA_H_VERT = 0, OMEGA_H_EDGE = 1, OMEGA_H_TRI = 2, OMEGA_H_TET = 3 };

enum Omega_h_Op { OMEGA_H_MIN, OMEGA_H_MAX, OMEGA_H_SUM };

enum Omega_h_Xfer {
  OMEGA_H_DONT_TRANSFER,
  OMEGA_H_INHERIT,
  OMEGA_H_LINEAR_INTERP,
  OMEGA_H_POINTWISE,
  OMEGA_H_CONSERVE,
  OMEGA_H_GLOBAL,
  OMEGA_H_LENGTH,
  OMEGA_H_QUALITY,
  OMEGA_H_METRIC,
  OMEGA_H_CONSERVE_R3D
};

enum { OMEGA_H_XFERS = OMEGA_H_CONSERVE_R3D + 1 };

enum Omega_h_Parting {
  OMEGA_H_ELEM_BASED,
  OMEGA_H_GHOSTED,
  OMEGA_H_VERT_BASED,
};

enum Omega_h_Comparison { OMEGA_H_SAME, OMEGA_H_MORE, OMEGA_H_DIFF };

enum Omega_h_Outflags {
  OMEGA_H_DONT_SAVE = 0x0,
  OMEGA_H_DO_SAVE   = 0x1,
  OMEGA_H_DONT_VIZ  = 0x0,
  OMEGA_H_DO_VIZ    = 0x2,
  OMEGA_H_SAVE_AND_VIZ = OMEGA_H_DO_SAVE | OMEGA_H_DO_VIZ
};

#ifdef __cplusplus
extern "C" {
#endif

void Omega_h_fail(char const* format, ...) __attribute__((noreturn));

void Omega_h_init_internal(int* argc, char*** argv, char const* head_desc);

inline static void Omega_h_init(int* argc, char*** argv) {
  Omega_h_init_internal(argc, argv, OMEGA_H_VERSION);
}

void Omega_h_finalize(void);

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
