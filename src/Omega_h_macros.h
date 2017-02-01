#ifndef OMEGA_H_MACROS_H
#define OMEGA_H_MACROS_H

#include <Omega_h_config.h>

#define OMEGA_H_STRINGIFY(s) #s
/* apparently you need two macros to make a string */
#define OMEGA_H_TOSTRING(s) OMEGA_H_STRINGIFY(s)

/* this block of preprocessor code creates a string
   literal describing the version and compile options
   used to configure this header.
   this is used to ensure that during user application
   compiles, the header included matches the library linked to.
   it is not foolproof, but better than nothing */
#define OMEGA_H_VERSION                                                        \
  OMEGA_H_TOSTRING(OMEGA_H_VERSION_MAJOR)                                      \
  "." OMEGA_H_TOSTRING(OMEGA_H_VERSION_MINOR) "." OMEGA_H_TOSTRING(            \
      OMEGA_H_VERSION_PATCH) "+" OMEGA_H_TOSTRING(defined(OMEGA_H_USE_MPI))    \
      OMEGA_H_TOSTRING(defined(OMEGA_H_USE_KOKKOS))                            \
          OMEGA_H_TOSTRING(defined(OMEGA_H_USE_OPENMP))                        \
              OMEGA_H_TOSTRING(defined(OMEGA_H_USE_CUDA))                      \
                  OMEGA_H_TOSTRING(defined(OMEGA_H_USE_ZLIB))                  \
                  OMEGA_H_TOSTRING(defined(OMEGA_H_USE_LIBMESHB))                  \
                  OMEGA_H_TOSTRING(defined(OMEGA_H_USE_EGADS))                  \
                      OMEGA_H_TOSTRING(defined(OMEGA_H_CHECK_BOUNDS))

#define OMEGA_H_PRAGMA(x) _Pragma(#x)

#if defined(__clang__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(clang system_header)
#elif defined(__GNUC__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(GCC system_header)
#endif

#endif
