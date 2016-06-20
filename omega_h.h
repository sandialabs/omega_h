#ifndef OMEGA_H_H
#define OMEGA_H_H

#include "omega_h_config.h"

/* MPI include block */

#ifdef OSH_USE_MPI

/* on BlueGene/Q the default install
 * defines MPICH2_CONST to empty if we
 * don't define it first, causing tons
 * of compile errors.
 *
 * in addition, the mpicxx.h header is full
 * of "const MPICH2_CONST", probably as a workaround
 * for the first mistake above.
 * as a result, properly defining MPICH2_CONST to const
 * also causes compile errors.
 * luckily, we can avoid including mpicxx.h with
 * MPICH_SKIP_MPICXX.
 */
#ifdef __bgq__
#define MPICH2_CONST const
#define MPICH_SKIP_MPICXX
#endif

/* have Clang diagnostics ignore everything
 * inside mpi.h
 */
#ifdef __clang__
#pragma clang system_header
#endif
#include <mpi.h>

#endif //OSH_USE_MPI

/* standard C headers required to parse this file */

#include <assert.h>

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
