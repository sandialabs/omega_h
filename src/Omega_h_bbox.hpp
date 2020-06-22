#ifndef OMEGA_H_BBOX_HPP
#define OMEGA_H_BBOX_HPP

#include <Omega_h_affine.hpp>
#include <Omega_h_vector.hpp>

namespace Omega_h {

class Mesh;

#ifdef OMEGA_H_USE_KOKKOS

template <Int dim>
struct BBox {
  OMEGA_H_INLINE BBox() {}
  OMEGA_H_INLINE BBox(Vector<dim> x) : min(x), max(x) {}
  OMEGA_H_INLINE BBox(Vector<dim> min_, Vector<dim> max_)
      : min(min_), max(max_) {}
  Vector<dim> min;
  Vector<dim> max;
  /* playing the volatile game again (see int128.hpp) */
  OMEGA_H_INLINE void operator=(BBox<dim> const& rhs) volatile {
    min = rhs.min;
    max = rhs.max;
  }
  OMEGA_H_INLINE BBox(BBox<dim> const& rhs) : min(rhs.min), max(rhs.max) {}
  OMEGA_H_INLINE BBox(const volatile BBox<dim>& rhs)
      : min(rhs.min), max(rhs.max) {}
};

#else

template <Int dim>
struct BBox {
  Vector<dim> min;
  Vector<dim> max;
  inline BBox() = default;
  OMEGA_H_INLINE BBox(Vector<dim> x) : min(x), max(x) {}
  OMEGA_H_INLINE BBox(Vector<dim> min_, Vector<dim> max_)
      : min(min_), max(max_) {}
  inline BBox(BBox const&) = default;
  inline BBox(BBox&&) = default;
  inline BBox& operator=(BBox const&) = default;
  inline BBox& operator=(BBox&&) = default;
};

#endif

template <Int dim>
OMEGA_H_INLINE BBox<dim> unite(BBox<dim> a, BBox<dim> b) {
  BBox<dim> c;
  for (Int i = 0; i < dim; ++i) {
    c.min[i] = min2(a.min[i], b.min[i]);
    c.max[i] = max2(a.max[i], b.max[i]);
  }
  return c;
}

template <Int dim>
OMEGA_H_INLINE bool are_close(BBox<dim> a, BBox<dim> b) {
  return are_close(a.min, b.min) && are_close(a.max, b.max);
}

template <Int dim>
OMEGA_H_INLINE BBox<dim> make_equilateral(BBox<dim> bbox) {
  auto maxl = reduce(bbox.max - bbox.min, maximum<Real>());
  for (Int i = 0; i < dim; ++i) bbox.max[i] = bbox.min[i] + maxl;
  return bbox;
}

template <Int dim>
OMEGA_H_INLINE Affine<dim> get_affine_from_bbox_into_unit(BBox<dim> bbox) {
  Vector<dim> s;
  for (Int i = 0; i < dim; ++i) s[i] = 1.0 / (bbox.max[i] - bbox.min[i]);
  Affine<dim> a;
  a.r = diagonal(s);
  a.t = -(a.r * bbox.min);
  return a;
}

template <Int dim>
BBox<dim> find_bounding_box(Reals coords);

template <Int dim>
BBox<dim> get_bounding_box(Mesh* mesh);

extern template BBox<1> find_bounding_box<1>(Reals coords);
extern template BBox<2> find_bounding_box<2>(Reals coords);
extern template BBox<3> find_bounding_box<3>(Reals coords);

extern template BBox<1> get_bounding_box<1>(Mesh* mesh);
extern template BBox<2> get_bounding_box<2>(Mesh* mesh);
extern template BBox<3> get_bounding_box<3>(Mesh* mesh);

}  // end namespace Omega_h

#endif
