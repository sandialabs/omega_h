#include "config.hpp"

#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <type_traits>
#include <chrono>
#include <memory>
#include <climits>
#include <vector>
#include <cfloat>

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
#include "metric.hpp"
#include "simplices.hpp"
#include "size.hpp"
#include "quality.hpp"
#include "loop.hpp"
#include "functors.hpp"
#include "array.hpp"
#include "access.hpp"
#include "repro.hpp"
#include "sort.hpp"
#include "scan.hpp"
#include "indset.hpp"
#include "map.hpp"
#include "align.hpp"
#include "adjacency.hpp"
#include "tag.hpp"
#include "mesh.hpp"
#include "hilbert.hpp"
