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

#ifdef USE_KOKKOS
#ifdef __clang__
#pragma clang system_header
#endif
#include <Kokkos_Core.hpp>
#endif

#ifdef USE_KOKKOS
#define INLINE KOKKOS_INLINE_FUNCTION
#define LAMBDA KOKKOS_LAMBDA
#else
#define INLINE inline
#define LAMBDA [=]
#endif

#define NORETURN(x) do { assert(0); return x; } while(0)

typedef std::uint8_t  U8;
typedef std::uint16_t U16;
typedef std::uint32_t U32;
typedef std::uint64_t U64;
typedef std::uintptr_t UInt;
typedef double Real;

#define EPSILON 1e-10
#define PI 3.141592653589793

#include "control.hpp"
#include "timer.hpp"
#include "algorithm.hpp"
#include "few.hpp"
#include "algebra.hpp"
#include "space.hpp"
#include "qr.hpp"
#include "loop.hpp"
#include "array.hpp"
#include "access.hpp"
#include "int128.hpp"
