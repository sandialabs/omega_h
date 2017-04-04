#ifndef HILBERT_HPP
#define HILBERT_HPP

#include "Omega_h_internal.hpp"

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
/* original code had unsigned int here, we are going to 64 bits */
typedef std::uint64_t coord_t;  // char,short,int for up to 8,16,32 bits per
                                // word
INLINE
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
    for (i = n - 1; i >= 0; i--)
      if (X[i] & Q)
        X[0] ^= P;  // invert
      else {
        t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
  }  // exchange
}
INLINE
void AxestoTranspose(coord_t* X, int b, int n)  // position, #bits, dimension
{
  coord_t M = coord_t(1) << (b - 1), P, Q, t;
  int i;
  // Inverse undo
  for (Q = M; Q > 1; Q >>= 1) {
    P = Q - 1;
    for (i = 0; i < n; i++)
      if (X[i] & Q)
        X[0] ^= P;  // invert
      else {
        t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
  }  // exchange
     // Gray encode
  for (i = 1; i < n; i++) X[i] ^= X[i - 1];
  t = 0;
  for (Q = M; Q > 1; Q >>= 1)
    if (X[n - 1] & Q) t ^= Q - 1;
  for (i = 0; i < n; i++) X[i] ^= t;
}

/* end verbatim code, what follows are omega_h helpers */

// un-transpose:
//                  in[0] = A D G J M
//                  in[1] = B E H K N
//                  in[2] = C F I L O
// becomes:
//                 out[0] = A B C D E
//                 out[1] = F G H I J
//                 out[2] = K L M N O
INLINE void untranspose(coord_t const in[], coord_t out[], int b, int n) {
  for (int i = 0; i < n; ++i) out[i] = 0;
  for (int i = 0; i < (b * n); ++i)
    out[i / b] |= (((in[i % n] >> (b - 1 - (i / n))) & 1) << (b - 1 - (i % b)));
}

/* output a permutation from sorted points to input
   points, such that their ordering reflects the
   traversal of a fine-scale Hilbert curve over
   the bounding box of the points */
LOs sort_coords(Reals coords, Int dim);

}  // end namespace hilbert

}  // end namespace Omega_h

#endif
