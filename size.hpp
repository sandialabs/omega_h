INLINE Real triangle_area(Few<Vector<2>, 3> p) {
  return cross(p[1] - p[0], p[2] - p[0]) / 2.0;
}

INLINE Real triangle_area(Few<Vector<3>, 3> p) {
  return norm(cross(p[1] - p[0], p[2] - p[0])) / 2.0;
}

INLINE Real tet_volume(Few<Vector<3>, 4> p) {
  return (cross(p[1] - p[0], p[2] - p[0]) * (p[3] - p[0])) / 6.0;
}

/* we are interested in formulas for area and volume that require
   only edge lengths as input, since the metric tensor
   is best at evaluating distance and angles,
   not directly transforming coordinates.
   some codes will store the metric's eigendecomposition,
   but then interpolation is ad-hoc.
   other codes store the metric tensor directly, but
   then do an eigendecomposition to get the transformation
   to metric space so they can transform the element first
   and then use the classic shape measure.
   we would like to avoid having to compute the eigendecomposition
   of the metric tensor at all, which may be possible if
   we use the metric directly for edge lengths and then
   use an edge-length-based shape measure...

   besides Wikipedia, this paper is a useful source:

   Kahan, William M.
   "What has the Volume of a Tetrahedron to do with Computer Programming Languages?"

   This paper not only gives formulas for triangle and
   tetrahedron size from edge lengths, but also shows
   that they are numerically accurate (at least backward stable).
   This is useful to us because in anisotropic meshes
   we will be stressing this code numerically.
*/

/* see the Kahan paper.
   this computes (u - v + w) in away that improves
   numerical conditioning over the naive expression */

INLINE Real facial_difference(Real u, Real v, Real w) {
  return (max2(u, w) - v) + min2(u, w);
}

INLINE Real triangle_area_kahan(Real u, Real v, Real w)
{
  return sqrt(
        (u + v + w) *
        facial_difference(v, w, u) *
        facial_difference(w, u, v) *
        facial_difference(u, v, w)
      ) / 4.0;
}

/*
edge pairs that are opposite one another in topology:
(u,U) (v,V) (w,W)
*/

INLINE Real tet_volume_kahan(
    Real U, Real V, Real W,
    Real u, Real v, Real w)
{
  Real X_U = facial_difference(w, U, v);
  Real Y_V = facial_difference(u, V, w);
  Real Z_W = facial_difference(v, W, u);
  Real X_V = facial_difference(U, v, w);
  Real Y_W = facial_difference(V, w, u);
  Real Z_U = facial_difference(W, u, v);
  Real X_W = facial_difference(v, w, U);
  Real Y_U = facial_difference(w, u, V);
  Real Z_V = facial_difference(u, v, W);
  Real X_0 = (U + v + w);
  Real Y_0 = (V + w + u);
  Real Z_0 = (W + u + v);
  Real X = X_U * X_0;
  Real Y = Y_V * Y_0;
  Real Z = Z_W * Z_0;
  Real x = X_V * X_W;
  Real y = Y_W * Y_U;
  Real z = Z_U * Z_V;
  Real xi     = sqrt(x * Y * Z);
  Real eta    = sqrt(y * Z * X);
  Real zeta   = sqrt(z * X * Y);
  Real lambda = sqrt(x * y * z);
  return sqrt((xi     + eta    + zeta   - lambda) *
              (lambda + xi     + eta    - zeta  ) *
              (eta    + zeta   + lambda - xi    ) *
              (zeta   + lambda + xi     - eta   ))
         / (192 * u * v * w);
}
