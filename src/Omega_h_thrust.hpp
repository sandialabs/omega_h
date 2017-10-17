#ifndef OMEGA_H_THRUST_HPP
#define OMEGA_H_THRUST_HPP

#include <Omega_h_macros.h>

OMEGA_H_SYSTEM_HEADER

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <thrust/device_ptr.h>
#include <thrust/sort.h>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#endif
