#ifndef OMEGA_H_HILBERT_HPP
#define OMEGA_H_HILBERT_HPP

#include <cstdint>

#include <Omega_h_affine.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_vector.hpp>

namespace Omega_h {

namespace hilbert {

/* The following code is verbatim plus one bug fix from the paper:

   Skilling, John. "Programming the Hilbert curve."
   BAYESIAN INFERENCE AND MAXIMUM ENTROPY METHODS IN SCIENCE AND ENGINEERING:
   23rd International Workshop on Bayesian Inference and
   Maximum Entropy Methods in Science and Engineering.
   Vol. 707. No. 1. AIP Publishing, 2004.
*/

//+++++++++++++++++++++++++++ PUBLIC-DOMAIN SOFTWARE ++++++++++++++++++++++++++
// Functions: TransposetoAxes AxestoTranspose
// Purpose:   Transform in-place between Hilbert transpose and geometrical axes
// Example:   b=5 bits for each of n=3 coordinates.
//            15-bit Hilbert integer = A B C D E F G H I J K L M N O  is stored
//            as its Transpose
//                   X[0] = A D G J M                X[2]|
//                   X[1] = B E H K N                    | /X[1]
//                   X[2] = C F I L O               axes |/
//                          high  low                    O------ X[0]
//            Axes are stored conventionally as b-bit integers
// Author:    John Skilling  20 Apr 2001 to 11 Oct 2003
//-----------------------------------------------------------------------------
/* Dan Ibanez: original code had unsigned int here, we are going to 64 bits */
typedef std::uint64_t coord_t;  // char,short,int for up to 8,16,32 bits per
                                // word
OMEGA_H_INLINE
void TransposetoAxes(coord_t* X, int b, int n)  // position, #bits, dimension
{
  coord_t N = coord_t(2) << (b - 1), P, Q, t;
  int i;
  // Gray code by H ^ (H/2)
  t = X[n - 1] >> 1;
  /* bug fix: original code has a (>= 0) here,
     but the loop body reads X[i-1]
                    V                       */
  for (i = n - 1; i > 0; i--) X[i] ^= X[i - 1];
  X[0] ^= t;
  // Undo excess work
  for (Q = 2; Q != N; Q <<= 1) {
    P = Q - 1;
    for (i = n - 1; i >= 0; i--) {
      if (X[i] & Q)
        X[0] ^= P;  // invert
      else {
        t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
    }
  }  // exchange
}

OMEGA_H_INLINE
void AxestoTranspose(coord_t* X, int b, int n)  // position, #bits, dimension
{
  coord_t M = coord_t(1) << (b - 1), P, Q, t;
  int i;
  // Inverse undo
  for (Q = M; Q > 1; Q >>= 1) {
    P = Q - 1;
    for (i = 0; i < n; i++) {
      if (X[i] & Q)
        X[0] ^= P;  // invert
      else {
        t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
    }
  }  // exchange
  // Gray encode
  for (i = 1; i < n; i++) X[i] ^= X[i - 1];
  t = 0;
  for (Q = M; Q > 1; Q >>= 1) {
    if (X[n - 1] & Q) t ^= Q - 1;
  }
  for (i = 0; i < n; i++) X[i] ^= t;
}

// un-transpose:
//                  in[0] = A D G J M
//                  in[1] = B E H K N
//                  in[2] = C F I L O
// becomes:
//                 out[0] = A B C D E
//                 out[1] = F G H I J
//                 out[2] = K L M N O
OMEGA_H_INLINE void untranspose(
    coord_t const in[], coord_t out[], int b, int n) {
  for (int i = 0; i < n; ++i) out[i] = 0;
  for (int i = 0; i < (b * n); ++i) {
    out[i / b] |= (((in[i % n] >> (b - 1 - (i / n))) & 1) << (b - 1 - (i % b)));
  }
}

/* Dan Ibanez: end verbatim code, what follows are omega_h helpers */

/* converts a floating-point spatial coordinate into an integral
   1D Hilbert coordinate on an implicit regular grid.
   The grid is defined by an affine transformation
   which maps real space vectors into a unit box,
   where the number of implicit grid cells along
   one axis of the unit box is (2^nbits)).

   After doing this, it converts the integral grid cell indices into
   a "one-dimensional" Hilbert space-filling-curve integer.
   In practice, this integer may be up to (64*dim) bits, so it is stored
   as several 64-bit integers.

   It is the user's responsibility to ensure that the affine transformation
   maps all possible input points such that the resulting point has all
   its coordinates in the range [0.0, 1.0].
   This function will clamp those coordinates for additional safety.
 */
template <Int dim>
OMEGA_H_INLINE Few<hilbert::coord_t, dim> from_spatial(
    Affine<dim> to_unit_box, Int nbits, Vector<dim> coord) {
  auto unit_box_coord = to_unit_box * coord;
  hilbert::coord_t X[dim];
  for (Int j = 0; j < dim; ++j) {
    /* this is more of an assert, and allows coordinates to be slightly
       outside the unit box without too severe consequences */
    auto zero_to_one_coord = clamp(unit_box_coord[j], 0.0, 1.0);
    auto zero_to_2eP_coord = zero_to_one_coord * std::exp2(Real(nbits));
    X[j] = hilbert::coord_t(zero_to_2eP_coord);
    /* some values will just graze the acceptable range
       (with proper floating point math they are exactly
        equal to 2^(nbits), and we'll be safe with (>=) in case
       floating point math is even worse than that. */
    if (X[j] >= (hilbert::coord_t(1) << nbits)) {
      X[j] = (hilbert::coord_t(1) << nbits) - 1;
    }
  }
  hilbert::AxestoTranspose(X, nbits, dim);
  Few<hilbert::coord_t, dim> Y;
  hilbert::untranspose(X, &Y[0], nbits, dim);
  return Y;
}

/* output a permutation from sorted points to input
   points, such that their ordering reflects the
   traversal of a fine-scale Hilbert curve over
   the bounding box of the points */
LOs sort_coords(Reals coords, Int dim);

}  // end namespace hilbert

}  // end namespace Omega_h

#endif
