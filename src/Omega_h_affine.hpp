#ifndef OMEGA_H_AFFINE_HPP
#define OMEGA_H_AFFINE_HPP

#include <Omega_h_matrix.hpp>

namespace Omega_h {

template <Int dim>
struct Affine {
  Tensor<dim> r;
  Vector<dim> t;
};

template <Int dim>
OMEGA_H_INLINE Vector<dim> operator*(Affine<dim> a, Vector<dim> v) {
  return (a.r * v) + a.t;
}

template <Int dim>
OMEGA_H_INLINE Affine<dim> invert(Affine<dim> a) {
  Affine<dim> ai;
  ai.r = invert(a.r);
  ai.t = -(ai.r * a.t);
  return ai;
}

}  // namespace Omega_h

#endif
