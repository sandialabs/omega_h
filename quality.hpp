/* There are several kinds of shape functions that have
   been used successfully for mesh adaptation.
   Each kind has a triangle and a tet variant.

   Symbols used below:
     Q_{tet} tetrahedral quality measure, [0,1]
     Q_{tri} triangular  quality measure, [0,1]
     V       tetrahedron volume
     A       triangle area
     A_i     area of triangle (i) of a tet
     l_i     length of edge (i) of a triangle or tet

   The first is the mean ratio measure, used by
   SCOREC and INRIA.
   For tets, the SCOREC implementation is the mean
   ratio cubed:

   Q_{tet} = 15552 V^2 / (\sum_{i=1}^6 l_i^2)^3

source:
   Li, Xiangrong, Mark S. Shephard, and Mark W. Beall.
   "3D anisotropic mesh adaptation by mesh modification."
   Computer methods in applied mechanics and engineering
   194.48 (2005): 4915-4950.

   The INRIA implementation for tets should be the
   actual mean ratio; the cube root of the SCOREC one:

   Q_{tet} = (36/(3^{1/3})) (V^{2/3}/(\sum_{i=1}^6 l_i^2))

   15552 = (36^3)/3

source:
   Loseille, Adrien, Victorien Menier, and Frederic Alauzet.
   "Parallel Generation of Large-size Adapted Meshes."
   Procedia Engineering 124 (2015): 57-69.
   (the normalization factor seems to have been misprinted
    here, it was not inverted like the rest of the formula)

   For triangles, the SCOREC variant is:

   Q_{tri} = 48 A^2 / (\sum_{i=1}^6 l_i^2)^2

   Another pair of measures that has been tried recently
   are some scale-invariant smooth measures associated with
   element stiffness matrix conditioning:

   Q_{tet} = V / (A_{rms})^{3/2}
   A_{rms} = \sqrt{(1/4)\sum_{i=1}^4 A_i^2}

   Q_{tri} = A / (l_{rms})^2
   l_{rms} = \sqrt{(1/3)\sum_{i=1}^3 l_i^2}
source:
   Shewchuk, J.
   "What is a good linear finite element?
    interpolation, conditioning, anisotropy,
    and quality measures"
   Proceedings of the 11th International Meshing Roundtable

   Some more treatment of quality measures can
   be found in
source:
   Liu, Anwei, and Barry Joe.
   "Relationship between tetrahedron shape measures."
   BIT Numerical Mathematics 34.2 (1994): 268-287.

   We will start off using the mean ratio cubed measures:
 */

INLINE Real triangle_mean_ratio_squared(Real a, Few<Real, 3> lsq) {
  Real s = 0.0;
  for (Int i = 0; i < 3; ++i)
    s += lsq[i];
  return 48 * square(a) / square(s);
}

INLINE Real tet_mean_ratio_cubed(Real v, Few<Real, 6> lsq) {
  Real s = 0.0;
  for (Int i = 0; i < 6; ++i)
    s += lsq[i];
  return 15552.0 * square(v) / cube(s);
}

INLINE Real real_triangle_quality(Few<Vector<2>, 3> p) {
  auto b = simplex_basis<2,2>(p);
  auto a = triangle_area(b);
  if (a < 0)
    return a;
  Few<Real, 3> lsq;
  lsq[0] = norm_squared(b[0]);
  lsq[1] = norm_squared(p[2] - p[1]);
  lsq[2] = norm_squared(b[1]);
  return sqrt(triangle_mean_ratio_squared(a, lsq));
}

INLINE Real real_tet_quality(Few<Vector<3>, 4> p) {
  auto b = simplex_basis<3,3>(p);
  auto v = tet_volume(b);
  if (v < 0)
    return v;
  Few<Real, 6> lsq;
  lsq[0] = norm_squared(b[0]);
  lsq[1] = norm_squared(p[2] - p[1]);
  lsq[2] = norm_squared(b[1]);
  lsq[3] = norm_squared(b[2]);
  lsq[4] = norm_squared(p[3] - p[1]);
  lsq[5] = norm_squared(p[3] - p[2]);
  return cbrt(tet_mean_ratio_cubed(v, lsq));
}

/* note that we will always use a constant metric tensor over the whole
   element to compute its quality, because that way we are computing
   the quality of the element after a single linear transformation which
   is guaranteed not to turn it inside out.
   this is why edge lengths are recomputed using the metric interpolated
   to the element centroid */

INLINE Real metric_triangle_quality(Few<Vector<2>, 3> p, Matrix<2,2> metric) {
  auto b = simplex_basis<2,2>(p);
  auto a = triangle_area(b);
  if (a < 0)
    return a;
  a *= sqrt(determinant(metric));
  Few<Real, 3> lsq;
  lsq[0] = metric_product(metric, b[0]);
  lsq[1] = metric_product(metric, p[2] - p[1]);
  lsq[2] = metric_product(metric, b[1]);
  return sqrt(triangle_mean_ratio_squared(a, lsq));
}

/* This paper:

   Loseille, Adrien, Victorien Menier, and Frederic Alauzet.
   "Parallel Generation of Large-size Adapted Meshes."
   Procedia Engineering 124 (2015): 57-69.

   Mentions using $\sqrt{\det(M)}$ to compute volume in metric space. */

INLINE Real metric_tet_quality(Few<Vector<3>, 4> p, Matrix<3,3> metric) {
  auto b = simplex_basis<3,3>(p);
  auto v = tet_volume(b);
  if (v < 0)
    return v;
  v *= sqrt(determinant(metric));
  Few<Real, 6> lsq;
  lsq[0] = metric_product(metric, b[0]);
  lsq[1] = metric_product(metric, p[2] - p[1]);
  lsq[2] = metric_product(metric, b[1]);
  lsq[3] = metric_product(metric, b[2]);
  lsq[4] = metric_product(metric, p[3] - p[1]);
  lsq[5] = metric_product(metric, p[3] - p[2]);
  return cbrt(tet_mean_ratio_cubed(v, lsq));
}

INLINE Real real_element_quality(Few<Vector<2>, 3> p) {
  return real_triangle_quality(p);
}

INLINE Real real_element_quality(Few<Vector<3>, 4> p) {
  return real_tet_quality(p);
}

INLINE Real metric_element_quality(Few<Vector<2>, 3> p, Matrix<2,2> metric) {
  return metric_triangle_quality(p, metric);
}

INLINE Real metric_element_quality(Few<Vector<3>, 4> p, Matrix<3,3> metric) {
  return metric_tet_quality(p, metric);
}

struct RealElementQualities {
  Reals coords;
  RealElementQualities(Mesh const* mesh):coords(mesh->coords()) {}
  template <Int neev>
  DEVICE Real measure(Few<LO, neev> v) const {
    auto p = gather_vectors<neev, neev - 1>(coords, v);
    return real_element_quality(p);
  }
};

struct MetricElementQualities {
  Reals coords;
  Reals metrics;
  MetricElementQualities(Mesh const* mesh):
    coords(mesh->coords()),
    metrics(mesh->get_array<Real>(VERT, "metric"))
  {}
  template <Int neev>
  DEVICE Real measure(Few<LO, neev> v) const {
    auto p = gather_vectors<neev, neev - 1>(coords, v);
    auto metric = average_metrics(gather_symms<neev,neev - 1>(metrics, v));
    return metric_element_quality(p, metric);
  }
};

Reals measure_qualities(Mesh* mesh, LOs a2e);
Reals measure_qualities(Mesh* mesh);
