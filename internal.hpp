#include "omega_h.hpp"

#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <type_traits>
#include <chrono>
#include <memory>
#include <climits>
#include <vector>
#include <cfloat>

#ifdef OSH_USE_ZLIB
#include <zlib.h>
#endif

/* in some cases involving __constant__ data,
 * we can't have an inline function marked as
 * __host__ at all (which KOKKOS_INLINE_FUNCTION does),
 * so we need a separate macro for things that
 * execute only on the device */

#ifdef OSH_USE_CUDA
#define DEVICE __device__ inline
#else
#define DEVICE inline
#endif

#define NORETURN(x) do { assert(0); return x; } while(0)

#ifdef OSH_USE_CUDA
#define CONSTANT __constant__
#else
#define CONSTANT
#endif

typedef std::int8_t  I8;
typedef std::int16_t I16;
typedef std::int32_t I32;
typedef std::int64_t I64;
typedef I32          Int;
typedef I32          LO;
typedef I64          GO;
typedef double       Real;

#define EPSILON 1e-10
#define PI 3.141592653589793
#define MANTISSA_BITS 52

#include "control.hpp"
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
#include "mesh.hpp"
#include "indset.hpp"
#include "metric.hpp"
#include "size.hpp"
#include "quality.hpp"
#include "construct.hpp"
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
#include "gmsh.hpp"
#include "repro.hpp"
#include "dist.hpp"
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
