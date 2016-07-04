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

namespace osh {

#include "protect.hpp"
#include "timer.hpp"
#include "algorithm.hpp"
#include "few.hpp"
#include "int128.hpp"
#include "traits.hpp"
#include "atomics.hpp"
#include "algebra.hpp"
#include "qr.hpp"
#include "space.hpp"
#include "polynomial.hpp"
#include "eigen.hpp"
#include "simplices.hpp"
#include "loop.hpp"
#include "functors.hpp"
#include "array.hpp"
#include "access.hpp"
#include "sort.hpp"
#include "scan.hpp"
#include "map.hpp"
#include "align.hpp"
#include "graph.hpp"
#include "adjacency.hpp"
#include "tag.hpp"
#include "comm.hpp"
#include "remotes.hpp"
#include "indset.hpp"
#include "metric.hpp"
#include "size.hpp"
#include "quality.hpp"
#include "vtk.hpp"
#include "bbox.hpp"
#include "hilbert.hpp"
#include "file.hpp"
#include "base64.hpp"
#include "simplify.hpp"
#include "box.hpp"
#include "surface.hpp"
#include "mark.hpp"
#include "classify.hpp"
#include "reorder.hpp"
#include "linpart.hpp"
#include "owners.hpp"
#include "migrate.hpp"
#include "bcast.hpp"
#include "unmap_mesh.hpp"
#include "ghost.hpp"
#include "inertia.hpp"
#include "bipart.hpp"
#include "refine_qualities.hpp"
#include "refine_topology.hpp"
#include "modify.hpp"
#include "refine.hpp"
#include "stacktrace.hpp"
#include "transfer.hpp"
#include "collapse.hpp"
#include "coarsen.hpp"
#include "laplace.hpp"
#include "adapt.hpp"
#include "swap2d.hpp"
#include "swap3d_tables.hpp"
#include "swap3d_loop.hpp"
#include "swap3d_choice.hpp"
#include "swap3d.hpp"
#include "swap.hpp"
#include "consistent.hpp"
#include "construct.hpp"
#include "xml.hpp"

}  // end namespace osh
