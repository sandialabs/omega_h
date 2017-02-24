#ifndef OMEGA_H_MACROS_H
#define OMEGA_H_MACROS_H

#include <Omega_h_config.h>

#define OMEGA_H_STRINGIFY(s) #s
/* apparently you need two macros to make a string */
#define OMEGA_H_TOSTRING(s) OMEGA_H_STRINGIFY(s)

#ifdef OMEGA_H_USE_MPI
#define OMEGA_H_VERSION_0 "1"
#else
#define OMEGA_H_VERSION_0 "0"
#endif

#ifdef OMEGA_H_USE_KOKKOS
#define OMEGA_H_VERSION_1 "1"
#else
#define OMEGA_H_VERSION_1 "0"
#endif

#ifdef OMEGA_H_USE_OPENMP
#define OMEGA_H_VERSION_2 "1"
#else
#define OMEGA_H_VERSION_2 "0"
#endif

#ifdef OMEGA_H_USE_CUDA
#define OMEGA_H_VERSION_3 "1"
#else
#define OMEGA_H_VERSION_3 "0"
#endif

#ifdef OMEGA_H_USE_ZLIB
#define OMEGA_H_VERSION_4 "1"
#else
#define OMEGA_H_VERSION_4 "0"
#endif

#ifdef OMEGA_H_USE_LIBMESHB
#define OMEGA_H_VERSION_5 "1"
#else
#define OMEGA_H_VERSION_5 "0"
#endif

#ifdef OMEGA_H_USE_EGADS
#define OMEGA_H_VERSION_6 "1"
#else
#define OMEGA_H_VERSION_6 "0"
#endif

#ifdef OMEGA_H_USE_SEACASEXODUS
#define OMEGA_H_VERSION_7 "1"
#else
#define OMEGA_H_VERSION_7 "0"
#endif

#ifdef OMEGA_H_CHECK_BOUNDS
#define OMEGA_H_VERSION_8 "1"
#else
#define OMEGA_H_VERSION_8 "0"
#endif

/* this block of preprocessor code creates a string
   literal describing the version and compile options
   used to configure this header.
   this is used to ensure that during user application
   compiles, the header included matches the library linked to.
   it is not foolproof, but better than nothing */
#define OMEGA_H_VERSION                                                        \
  OMEGA_H_TOSTRING(OMEGA_H_VERSION_MAJOR)                                      \
  "." OMEGA_H_TOSTRING(OMEGA_H_VERSION_MINOR) "." OMEGA_H_TOSTRING(            \
      OMEGA_H_VERSION_PATCH) "+" OMEGA_H_VERSION_0  OMEGA_H_VERSION_1                            \
      OMEGA_H_VERSION_2 OMEGA_H_VERSION_3 OMEGA_H_VERSION_4                  \
      OMEGA_H_VERSION_5 OMEGA_H_VERSION_6 OMEGA_H_VERSION_7                  \
      OMEGA_H_VERSION_8

#define OMEGA_H_PRAGMA(x) _Pragma(#x)

#if defined(__clang__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(clang system_header)
#elif defined(__GNUC__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(GCC system_header)
#endif

#endif
