template <Int dim>
Read<I64> hilbert_dists_from_coords(Reals coords) {
  auto bbox = find_bounding_box<dim>(coords);
  Real maxl = 0;
  for (Int i = 0; i < dim; ++i)
    maxl = max2(maxl, bbox.max[i] - bbox.min[i]);
  int expo;
  frexp(maxl, &expo);
  Real unit = exp2(Real(expo - MANTISSA_BITS));
  LO npts = coords.size() / dim;
  Write<I64> out;
  auto f = LAMBDA(LO i) {
    hilbert::coord_t X[dim];
    /* floating-point coordinate to fine-grid integer coordinate,
       should be non-negative since we subtract the BBox min */
    for (I8 j = 0; j < dim; ++j) {
      X[j] = hilbert::coord_t((coords[i * dim + j] - bbox.min[j]) / unit);
      CHECK(X[j] < (hilbert::coord_t(1) << MANTISSA_BITS));
    }
    hilbert::AxestoTranspose(X, MANTISSA_BITS, dim);
    hilbert::coord_t Y[dim];
    hilbert::untranspose(X, Y, MANTISSA_BITS, dim);
    for (I8 j = 0; j < dim; ++j)
    /* this cast *should* be safe... */
      out[i * dim + j] = static_cast<I64>(Y[j]);
  };
  parallel_for(npts, f);
  return out;
}

template Read<I64> hilbert_dists_from_coords<2>(Reals coords);
template Read<I64> hilbert_dists_from_coords<3>(Reals coords);
