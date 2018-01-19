#ifndef OMEGA_H_DEFINES_HPP
#define OMEGA_H_DEFINES_HPP

#include <cstdint>

#include <Omega_h_c.h>

namespace Omega_h {

typedef std::int8_t I8;
typedef std::int16_t I16;
typedef std::int32_t I32;
typedef std::int64_t I64;
typedef I8 Byte;
typedef I32 Int;
typedef I32 LO;
typedef I32 ClassId;
typedef I64 GO;
typedef double Real;

constexpr Real PI = OMEGA_H_PI;
constexpr Real EPSILON = OMEGA_H_EPSILON;
constexpr Int MANTISSA_BITS = OMEGA_H_MANTISSA_BITS;

enum { DIMS = OMEGA_H_DIMS };

enum : Int {
  VERT = OMEGA_H_VERT,
  EDGE = OMEGA_H_EDGE,
  FACE = OMEGA_H_FACE,
  REGION = OMEGA_H_REGION
};

}  // end namespace Omega_h

#endif
