#ifndef BBOX_HPP
#define BBOX_HPP

#include "algebra.hpp"

namespace Omega_h {

template <Int dim>
INLINE BBox<dim> unite(BBox<dim> a, BBox<dim> b) {
  BBox<dim> c;
  for (Int i = 0; i < dim; ++i) {
    c.min[i] = min2(a.min[i], b.min[i]);
    c.max[i] = max2(a.max[i], b.max[i]);
  }
  return c;
}

template <Int dim>
INLINE bool are_close(BBox<dim> a, BBox<dim> b) {
  return are_close(a.min, b.min) && are_close(a.max, b.max);
}

template <Int dim>
BBox<dim> find_bounding_box(Reals coords);

extern template BBox<2> find_bounding_box<2>(Reals coords);
extern template BBox<3> find_bounding_box<3>(Reals coords);

}  // end namespace Omega_h

#endif
