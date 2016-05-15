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

#include "kokkos.hpp"

#define NORETURN(x) do { assert(0); return x; } while(0)

typedef std::uint8_t  U8;
typedef std::uint16_t U16;
typedef std::uint32_t U32;
typedef std::uint64_t U64;
typedef U64           UInt;
typedef U64           GO;
typedef U32           LO;
typedef double Real;

#define EPSILON 1e-10
#define PI 3.141592653589793
#define MANTISSA_BITS 52

#include "control.hpp"
#include "timer.hpp"
#include "algorithm.hpp"
#include "few.hpp"
#include "int128.hpp"
#include "traits.hpp"
#include "algebra.hpp"
#include "qr.hpp"
#include "space.hpp"
#include "polynomial.hpp"
#include "size.hpp"
#include "quality.hpp"
#include "loop.hpp"
#include "functors.hpp"
#include "array.hpp"
#include "access.hpp"
#include "repro.hpp"
#include "sort.hpp"
