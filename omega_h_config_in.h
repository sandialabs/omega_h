#ifndef OMEGA_H_HPP
#define OMEGA_H_HPP

#define OSH_MAJOR 0
#define OSH_MINOR 5

#cmakedefine OSH_USE_MPI
#cmakedefine OSH_USE_KOKKOS
#cmakedefine OSH_USE_OPENMP
#cmakedefine OSH_USE_CUDA
#cmakedefine OSH_USE_ZLIB
#cmakedefine OSH_CHECK_BOUNDS

/* this block of preprocessor code creates a string
   literal describing the version and compile options
   used to configure this header.
   this is used to ensure that during user application
   compiles, the header included matches the library linked to.
   it is not foolproof, but better than nothing */
#define OSH_TOSTR2(s) #s
#define OSH_TOSTR(s) OSH_TOSTR2(s)
#ifdef OSH_USE_MPI
#define OSH_MPI_STR "MPI"
#else
#define OSH_MPI_STR ""
#endif
#ifdef OSH_USE_KOKKOS
#define OSH_KOKKOS_STR "Kokkos"
#else
#define OSH_KOKKOS_STR ""
#endif
#ifdef OSH_USE_OPENMP
#define OSH_OPENMP_STR "OpenMP"
#else
#define OSH_OPENMP_STR ""
#endif
#ifdef OSH_USE_CUDA
#define OSH_CUDA_STR "CUDA"
#else
#define OSH_CUDA_STR ""
#endif
#ifdef OSH_USE_ZLIB
#define OSH_ZLIB_STR "zlib"
#else
#define OSH_ZLIB_STR ""
#endif
#define OSH_DESC \
"omega_h v" OSH_TOSTR(OSH_MAJOR) "." OSH_TOSTR(OSH_MINOR) \
": " OSH_MPI_STR "," OSH_KOKKOS_STR "," \
OSH_OPENMP_STR "," OSH_CUDA_STR "," OSH_ZLIB_STR

#endif
