#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include "algebra.hpp"

namespace Omega_h {

template <Int degree>
struct Roots {
  Int n;  // number of unique roots
  Few<Real, degree> values;
  Few<Int, degree> mults;  // multiplicities
};

// solve cubic equation x^3 + a_2 * x^2 + a_1 * x + a_0 = 0
// this code assumes that the solution does not have complex roots !
// the return value is the number of distinct roots,
// the two output arrays contain root values and multiplicities.
// roots within an *absolute* distance of (eps) are considered
// the same.
INLINE Roots<3> find_polynomial_roots(Few<Real, 3> coeffs, Real eps = 1e-4) {
  auto a_0 = coeffs[0];
  auto a_1 = coeffs[1];
  auto a_2 = coeffs[2];
  Few<Real, 3> roots;
  Few<Int, 3> mults;
  // http://mathworld.wolfram.com/CubicFormula.html
  Real p = (3. * a_1 - square(a_2)) / 3.;
  Real q = (9. * a_1 * a_2 - 27. * a_0 - 2. * cube(a_2)) / 27.;
  Real Q = p / 3.;
  Real R = q / 2.;
  Real D = cube(Q) + square(R);
  Real shift = -a_2 / 3.;
  if (D >= 0.0) {
    Real S = cbrt(R + sqrt(D));
    Real T = cbrt(R - sqrt(D));
    Real B = S + T;
    Real z_1 = shift + B;
    Real z_23_real = shift - (1. / 2.) * B;
    // warning: we simply assume the imaginary component is small !
    // Real z_23_imag = (1. / 2.) * sqrt(3.) * A;
    roots[0] = z_1;
    roots[1] = roots[2] = z_23_real;
  } else {
    // D < 0 implies Q < 0, since R^2 must be positive
    auto cos_theta = R / sqrt(-cube(Q));
    Real theta = acos(clamp(cos_theta, 1.0, -1.0));
    Real radius = 2. * sqrt(-Q);
    Real z_1 = radius * cos((theta) / 3.) + shift;
    Real z_2 = radius * cos((theta + 2. * PI) / 3.) + shift;
    Real z_3 = radius * cos((theta - 2. * PI) / 3.) + shift;
    roots[0] = z_1;
    roots[1] = z_2;
    roots[2] = z_3;
  }
  mults[0] = mults[1] = mults[2] = 1;
  // post-processing, decide whether pairs of roots are
  // close enough to be called repeated roots
  // first step, if two roots are close, then
  // move them to the second and third slots
  if (fabs(roots[0] - roots[1]) < eps) {
    swap2(roots[0], roots[2]);
  } else if (fabs(roots[0] - roots[2]) < eps) {
    swap2(roots[0], roots[1]);
  } else if (fabs(roots[1] - roots[2]) < eps) {
    // no need to swap, they're already there
  } else {
    // no pairs were close, all three roots are distinct
    return {3, roots, mults};
  }
  // if we're here, two close roots are in [1] and [2]
  roots[1] = average(roots[1], roots[2]);
  mults[1] = 2;
  // lets see if they are all the same
  if (fabs(roots[0] - roots[1]) < eps) {
    // roots[1] is already an average, weight it properly
    roots[0] = (1. / 3.) * roots[0] + (2. / 3.) * roots[1];
    mults[0] = 3;
    return {1, roots, mults};
  }
  return {2, roots, mults};
}

// solve quadratic equation x^2 + a * x + b = 0
INLINE Roots<2> find_polynomial_roots(Few<Real, 2> coeffs, Real eps = 1e-4) {
  auto a = coeffs[1];
  auto b = coeffs[0];
  Few<Real, 2> roots;
  Few<Int, 2> mults;
  Real disc = square(a) - 4. * b;
  if (fabs(disc) < eps) {
    mults[0] = 2;
    roots[0] = -a / 2.;
    roots[1] = roots[0];
    mults[1] = mults[0];
    return {1, roots, mults};
  }
  if (disc > 0.0) {
    mults[0] = 1;
    mults[1] = 1;
    roots[0] = (-a + sqrt(disc)) / 2.;
    roots[1] = (-a - sqrt(disc)) / 2.;
    return {2, roots, mults};
  }
  return {0, roots, mults};
}

}  // end namespace Omega_h

#endif
