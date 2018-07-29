#include "Omega_h_bbox.hpp"

#include "Omega_h_reduce.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

template <Int dim>
struct BBoxFunctor {
  using value_type = BBox<dim>;
  Reals coords_;
  BBoxFunctor(Reals coords) : coords_(coords) {}
  OMEGA_H_INLINE void init(value_type& update) const {
    for (Int i = 0; i < dim; ++i) {
      update.min[i] = ArithTraits<Real>::max();
      update.max[i] = ArithTraits<Real>::min();
    }
  }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = unite(update, input);
  }
  OMEGA_H_DEVICE void operator()(Int i, value_type& update) const {
    update = unite(update, BBox<dim>(get_vector<dim>(coords_, i)));
  }
};

template <Int dim>
BBox<dim> find_bounding_box(Reals coords) {
  auto npts = divide_no_remainder(coords.size(), dim);
  return parallel_reduce(npts, BBoxFunctor<dim>(coords), "find_bounding_box");
}

template BBox<1> find_bounding_box<1>(Reals coords);
template BBox<2> find_bounding_box<2>(Reals coords);
template BBox<3> find_bounding_box<3>(Reals coords);

template <Int dim>
BBox<dim> get_bounding_box(Mesh* mesh) {
  auto bb = find_bounding_box<dim>(mesh->coords());
  for (Int i = 0; i < dim; ++i) {
    bb.min[i] = mesh->comm()->allreduce(bb.min[i], OMEGA_H_MIN);
    bb.max[i] = mesh->comm()->allreduce(bb.max[i], OMEGA_H_MAX);
  }
  return bb;
}

template BBox<1> get_bounding_box<1>(Mesh* mesh);
template BBox<2> get_bounding_box<2>(Mesh* mesh);
template BBox<3> get_bounding_box<3>(Mesh* mesh);

}  // end namespace Omega_h
