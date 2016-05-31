#include "config.hpp"

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

#ifdef USE_ZLIB
#include <zlib.h>
#endif

#include "kokkos.hpp"

#define NORETURN(x) do { assert(0); return x; } while(0)

#ifdef USE_CUDA
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
#include "indset.hpp"
#include "map.hpp"
#include "align.hpp"
#include "graph.hpp"
#include "adjacency.hpp"
#include "tag.hpp"
#include "comm.hpp"
#include "remotes.hpp"
#include "mesh.hpp"
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
