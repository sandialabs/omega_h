#include "internal.hpp"

namespace osh {

namespace hilbert {

/* for each set of (dim) floating-point coordinates, this function
   outputs a set of (dim) 64-bit integers which represent the
   closest point of a fine-grid Hilbert curve to the coordinates.
   the resolution of the grid is chosen to be 52 bits (the floating-point
   mantissa size), giving 2^52 grid points per axis,
   and is scaled to the bounding box of the coordinates.
   the output integers are such that sort_by_keys() will sort along the
   Hilbert curve.
   More precisely, the bits of the Hilbert distance are spread evenly
   among the integers, with the first integer getting the most significant
   bits, and the last integer getting the least significant bits. */

template <Int dim>
static Read<I64> dists_from_coords_dim(Reals coords) {
  auto bbox = find_bounding_box<dim>(coords);
  Real maxl = 0;
  for (Int i = 0; i < dim; ++i) maxl = max2(maxl, bbox.max[i] - bbox.min[i]);
  LO npts = coords.size() / dim;
  Write<I64> out(npts * dim);
  auto f = LAMBDA(LO i) {
    hilbert::coord_t X[dim];
    Int nbits = MANTISSA_BITS;
    for (Int j = 0; j < dim; ++j) {
      /* floating-point coordinate to fine-grid integer coordinate,
         should be non-negative since we subtract the BBox min */
      Real coord = coords[i * dim + j];
      Real zero_to_one_coord = (coord - bbox.min[j]) / maxl;
      Real zero_to_2eP_coord = zero_to_one_coord * exp2(Real(nbits));
      X[j] = hilbert::coord_t(zero_to_2eP_coord);
      /* some values will just graze the acceptable range
         (with proper floating point math they are exactly
          equal to 2^(nbits), and we'll be safe with (>=) in case
         floating point math is even worse than that. */
      if (X[j] >= (hilbert::coord_t(1) << nbits))
        X[j] = (hilbert::coord_t(1) << nbits) - 1;
    }
    hilbert::AxestoTranspose(X, nbits, dim);
    hilbert::coord_t Y[dim];
    hilbert::untranspose(X, Y, nbits, dim);
    for (Int j = 0; j < dim; ++j) /* this cast *should* be safe... */
      out[i * dim + j] = static_cast<I64>(Y[j]);
  };
  parallel_for(npts, f);
  return out;
}

static Read<I64> dists_from_coords(Reals coords, Int dim) {
  if (dim == 3) return dists_from_coords_dim<3>(coords);
  if (dim == 2) return dists_from_coords_dim<2>(coords);
  NORETURN(Read<I64>());
}

LOs sort_coords(Reals coords, Int dim) {
  auto keys = hilbert::dists_from_coords(coords, dim);
  return sort_by_keys(keys, dim);
}

}  // end namespace hilbert

} //end namespace osh
