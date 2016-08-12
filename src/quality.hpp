#ifndef QUALITY_HPP
#define QUALITY_HPP

#include "algebra.hpp"
#include "host_few.hpp"
#include "size.hpp"

namespace osh {

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

INLINE Real perfect_tri_area() {
  return sqrt(3.) / 4.;
}

INLINE Real tri_mean_ratio(Real a, Real msl) {
  return a / (perfect_tri_area() * msl);
}

INLINE Real perfect_tet_volume_squared() {
  return 1. / 72.;
}

INLINE Real tet_mean_ratio_cubed(Real v, Real msl) {
  return square(v) / (perfect_tet_volume_squared() * cube(msl));
}

INLINE Real tet_mean_ratio(Real v, Real msl) {
  return cbrt(tet_mean_ratio_cubed(v, msl));
}

template <Int dim>
struct MeanRatio;

template <>
struct MeanRatio<2> {
  INLINE static Real get(Real a, Real msl) {
    return tri_mean_ratio(a, msl);
  }
};

template <>
struct MeanRatio<3> {
  INLINE static Real get(Real v, Real msl) {
    return tet_mean_ratio(v, msl);
  }
};

template <Int dim>
INLINE Real mean_ratio(Real v, Real msl) {
  return MeanRatio<dim>::get(v, msl);
}

template <typename EdgeVectors>
Real mean_squared_real_length(EdgeVectors edge_vectors) {
  auto nedges = EdgeVectors::size;
  Real msl = 0;
  for (Int i = 0; i < nedges; ++i) {
    msl += norm_squared(edge_vectors[i]);
  }
  return msl / nedges;
}

template <Int dim>
INLINE Real real_element_quality(Few<Vector<dim>, dim + 1> p) {
  auto b = simplex_basis<dim, dim>(p);
  auto s = element_size(b);
  if (s < 0) return s;
  auto ev = element_edge_vectors(p, b);
  auto msl = mean_squared_real_length(ev);
  return mean_ratio<dim>(s, msl);
}

/* note that we will always use a constant metric tensor over the whole
   element to compute its quality, because that way we are computing
   the quality of the element after a single linear transformation which
   is guaranteed not to turn it inside out.
   this is why edge lengths are recomputed using the metric interpolated
   to the element centroid */

/* This paper:

   Loseille, Adrien, Victorien Menier, and Frederic Alauzet.
   "Parallel Generation of Large-size Adapted Meshes."
   Procedia Engineering 124 (2015): 57-69.

   Mentions using $\sqrt{\det(M)}$ to compute volume in metric space. */

template <Int dim, typename EdgeVectors>
Real mean_squared_metric_length(EdgeVectors edge_vectors, Matrix<dim, dim> metric) {
  auto nedges = EdgeVectors::size;
  Real msl = 0;
  for (Int i = 0; i < nedges; ++i) {
    msl += metric_product(metric, edge_vectors[i]);
  }
  return msl / nedges;
}

template <Int dim>
INLINE Real metric_element_quality(Few<Vector<dim>, dim + 1> p,
                                   Matrix<dim, dim> metric) {
  auto b = simplex_basis<dim, dim>(p);
  auto s = element_size(b) * sqrt(determinant(metric));
  if (s < 0) return s;
  auto ev = element_edge_vectors(p, b);
  auto msl = mean_squared_metric_length(ev, metric);
  return mean_ratio<dim>(s, msl);
}

struct RealElementQualities {
  Reals coords;
  RealElementQualities(Mesh const* mesh) : coords(mesh->coords()) {}
  template <Int neev>
  DEVICE Real measure(Few<LO, neev> v) const {
    auto p = gather_vectors<neev, neev - 1>(coords, v);
    return real_element_quality(p);
  }
};

struct MetricElementQualities {
  Reals coords;
  Reals metrics;
  MetricElementQualities(Mesh const* mesh)
      : coords(mesh->coords()),
        metrics(mesh->get_array<Real>(VERT, "metric")) {}
  template <Int neev>
  DEVICE Real measure(Few<LO, neev> v) const {
    auto p = gather_vectors<neev, neev - 1>(coords, v);
    auto metric = average_metrics(gather_symms<neev, neev - 1>(metrics, v));
    return metric_element_quality(p, metric);
  }
};

Reals measure_qualities(Mesh* mesh, LOs a2e);
Reals measure_qualities(Mesh* mesh);

constexpr Int nquality_histogram_buckets = 10;

typedef HostFew<GO, nquality_histogram_buckets> QualityHistogram;

QualityHistogram get_quality_histogram(Mesh* mesh);
void print_quality_histogram(QualityHistogram histogram);

}  // end namespace osh

#endif
