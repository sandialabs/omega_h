#ifndef OMEGA_H_MACROS_H
#define OMEGA_H_MACROS_H

#include <Omega_h_config.h>

#define OMEGA_H_STRINGIFY(s) #s
/* apparently you need two macros to make a string */
#define OMEGA_H_TOSTRING(s) OMEGA_H_STRINGIFY(s)

#define OMEGA_H_PRAGMA(x) _Pragma(#x)

#if defined(__clang__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(clang system_header)
#elif defined(__GNUC__)
#define OMEGA_H_SYSTEM_HEADER OMEGA_H_PRAGMA(GCC system_header)
#endif

#define OMEGA_H_PI 3.14159265358979323
//                  ^    ^    ^    ^
//                  0    5   10   15
#define OMEGA_H_EPSILON 1e-10
#define OMEGA_H_MANTISSA_BITS 52

#if defined(__clang__)
#define OMEGA_H_FALLTHROUGH [[clang::fallthrough]]
#else
#define OMEGA_H_FALLTHROUGH ((void)0)
#endif

#endif
