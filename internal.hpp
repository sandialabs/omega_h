#ifndef INTERNAL_HPP
#define INTERNAL_HPP

#include "omega_h.hpp"

/* C++ standard includes */

#include <array>
#include <cerrno>
#include <cfloat>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>
#include <type_traits>
#include <utility>
#include <map>

/* C++ ABI and POSIX includes */

#include <cxxabi.h>
#include <execinfo.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>

/* Third party libraries */

#ifdef OSH_USE_ZLIB
#include <zlib.h>
#endif

#if defined(OSH_USE_CUDA)
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#elif defined(OSH_USE_OPENMP)
#include <algorithm>
#include <omp.h>
#include "intel_sort/pss_common.hpp"
#include "intel_sort/parallel_stable_sort.hpp"
#else
#include <algorithm>
#endif

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
#define PI 3.141592653589793
#define MANTISSA_BITS 52

#endif
