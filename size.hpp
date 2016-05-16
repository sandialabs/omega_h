template <Int sdim, Int edim>
INLINE Few<Vector<sdim>, edim> simplex_basis(Few<Vector<sdim>, edim + 1> p) {
  Few<Vector<sdim>, edim> b;
  for (Int i = 0; i < edim; ++i)
    b[i] = p[i + 1] - p[0];
  return b;
}

INLINE Real triangle_area(Few<Vector<2>, 2> b) {
  return cross(b[0], b[1]) / 2.0;
}

INLINE Real triangle_area(Few<Vector<3>, 2> b) {
  return norm(cross(b[0], b[1])) / 2.0;
}

INLINE Real tet_volume(Few<Vector<3>, 3> b) {
  return (cross(b[0], b[1]) * b[2]) / 6.0;
}

/* the following edge length functions could be replaced
   with loops and lookup tables, but last time I did
   that things got noticeably slower */
template <Int dim>
INLINE Few<Real, 3> triangle_edge_lengths_squared(
    Few<Vector<dim>, 3> p,
    Few<Vector<dim>, 2> b) {
  Few<Real, 3> lsq;
  lsq[0] = norm_squared(b[0]);
  lsq[1] = norm_squared(p[2] - p[1]);
  lsq[2] = norm_squared(b[1]);
  return lsq;
}

template <Int dim>
INLINE Few<Real, 6> tet_edge_lengths_squared(
    Few<Vector<dim>, 4> p,
    Few<Vector<dim>, 3> b) {
  Few<Real, 6> lsq;
  lsq[0] = norm_squared(b[0]);
  lsq[1] = norm_squared(p[2] - p[1]);
  lsq[2] = norm_squared(b[1]);
  lsq[3] = norm_squared(b[2]);
  lsq[4] = norm_squared(p[3] - p[1]);
  lsq[5] = norm_squared(p[3] - p[2]);
  return lsq;
}
