#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <type_traits>
#include <chrono>
#include <memory>
typedef std::uint8_t  U8;
typedef std::uint16_t U16;
typedef std::uint32_t U32;
typedef std::uint64_t U64;
typedef std::uintptr_t UInt;
typedef double Real;
#define INLINE inline
#define EPSILON 1e-10
#define PI 3.141592653589793
#define NORETURN(x) do { assert(0); return x; } while(0)
#define LAMBDA [=]
#include "control.hpp"
#include "timer.hpp"
#include "algorithm.hpp"
#include "few.hpp"
#include "algebra.hpp"
#include "space.hpp"
#include "qr.hpp"
#include "loop.hpp"
#include "array.hpp"
