#ifndef OMEGA_H_R3D_HPP
#define OMEGA_H_R3D_HPP

#include <Omega_h_vector.hpp>

#ifdef OMEGA_H_USE_CUDA
#define R3D_USE_CUDA
#endif

#include <r3d.hpp>

namespace Omega_h {

template <Int n>
OMEGA_H_INLINE Vector<n> from_r3d(r3d::Vector<n> a) {
  Vector<n> b;
  for (Int i = 0; i < n; ++i) b[i] = a[i];
  return b;
}

template <Int n, Int dim>
OMEGA_H_INLINE Few<Vector<dim>, n> from_r3d(r3d::Few<r3d::Vector<dim>, n> a) {
  Few<Vector<dim>, n> b;
  for (Int i = 0; i < n; ++i) b[i] = from_r3d(a[i]);
  return b;
}

template <Int n>
OMEGA_H_INLINE r3d::Vector<n> to_r3d(Vector<n> a) {
  r3d::Vector<n> b;
  for (Int i = 0; i < n; ++i) b[i] = a[i];
  return b;
}

template <Int n, Int dim>
OMEGA_H_INLINE r3d::Few<r3d::Vector<dim>, n> to_r3d(Few<Vector<dim>, n> a) {
  r3d::Few<r3d::Vector<dim>, n> b;
  for (Int i = 0; i < n; ++i) b[i] = to_r3d(a[i]);
  return b;
}

}  // end namespace Omega_h

#endif
