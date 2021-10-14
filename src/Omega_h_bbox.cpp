#include "Omega_h_bbox.hpp"

#include "Omega_h_int_iterator.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_reduce.hpp"

namespace Omega_h {

template <Int dim>
struct UniteOp {
  OMEGA_H_INLINE BBox<dim> operator()(BBox<dim> a, BBox<dim> b) const {
    return unite(a, b);
  }
};

template <Int dim>
struct GetBBoxOp {
  Reals coords;
  GetBBoxOp(Reals coords_in) : coords(coords_in) {}
  OMEGA_H_INLINE BBox<dim> operator()(LO i) const {
    return BBox<dim>(get_vector<dim>(coords, i));
  }
};




//find the bbox enclosing all listed points
template <Int dim>
BBox<dim> find_bounding_box(Reals coords) {
  auto npts = divide_no_remainder(coords.size(), dim);
  BBox<dim> init;
  for (Int i = 0; i < dim; ++i) {
    init.min[i] = ArithTraits<Real>::max();
    init.max[i] = ArithTraits<Real>::min();
  }
#if defined(OMEGA_H_USE_KOKKOS) and !defined(OMEGA_H_USE_CUDA) and !defined(OMEGA_H_USE_OPENMP)
  BBox<dim> res;
//  const auto transform = GetBBoxOp<dim>(coords);
//  using space = typename Kokkos::View<int*>::memory_space; //HACK
//  if (npts > 0) {
//    Kokkos::parallel_reduce(
//      Kokkos::RangePolicy<>(0, npts),
//      KOKKOS_LAMBDA(int i, BBox<dim>& update) {
//        update = transform(i);
//      }, BBoxUnion<space,dim>(res));
//  }
//  return res;
#else
  return transform_reduce(IntIterator(0), IntIterator(npts), init,
      UniteOp<dim>(), GetBBoxOp<dim>(coords));
#endif
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
