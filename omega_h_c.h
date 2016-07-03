#ifndef OMEGA_H_C_H
#define OMEGA_H_C_H

#include "omega_h_config.h"

#define OSH_PRAGMA(x) _Pragma(#x)

#if defined(__clang__)
#define OSH_SYSTEM_HEADER OSH_PRAGMA(clang system_header)
#elif defined(__GNUC__)
#define OSH_SYSTEM_HEADER OSH_PRAGMA(GCC system_header)
#endif

#include "omega_h_mpi.h"

/* standard C headers required to parse this file */

#include <assert.h>

enum osh_type {
  OSH_I8  = 0,
  OSH_I32 = 2,
  OSH_I64 = 3,
  OSH_F64 = 5,
};

enum osh_op {
  OSH_MIN,
  OSH_MAX,
  OSH_SUM
};

enum osh_xfer {
  OSH_DONT_TRANSFER,
  OSH_INHERIT,
  OSH_LINEAR_INTERP,
  OSH_POINTWISE,
  OSH_CONSERVE,
  OSH_GLOBAL,
  OSH_LENGTH,
  OSH_QUALITY,
  OSH_METRIC
};

enum osh_parting {
  OSH_ELEM_BASED,
  OSH_GHOSTED,
  OSH_VERT_BASED,
};

enum osh_comparison {
  OSH_SAME,
  OSH_MORE,
  OSH_DIFF
};

#ifdef __cplusplus
extern "C" {
#endif

void osh_fail(char const* format, ...) __attribute__((noreturn));

void osh_init_internal(int* argc, char*** argv, char const* head_desc);

inline static void osh_init(int* argc, char*** argv) {
  osh_init_internal(argc, argv, OSH_VERSION);
}

void osh_finalize(void);

#ifdef __cplusplus
} //end of extern "C" block
#endif

#ifdef __CUDA_ARCH__
#define OSH_CHECK(cond) assert(cond)
#else
#define OSH_CHECK(cond) ((cond) ? ((void)0) : \
  osh_fail("assertion %s failed at %s +%d\n", #cond, __FILE__, __LINE__))
#endif

#endif
