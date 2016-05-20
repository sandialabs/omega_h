template <Int dim>
struct BBoxFunctor {
  typedef BBox<dim> value_type;
  Reals coords_;
  BBoxFunctor(Reals coords):coords_(coords) {}
  INLINE void init(value_type& update) const {
    for (Int i = 0; i < dim; ++i) {
      update.min[i] = ArithTraits<Real>::max();
      update.max[i] = ArithTraits<Real>::min();
    }
  }
  INLINE void join(volatile value_type& update,
      const volatile value_type& input) const {
    update = unite(update, input);
  }
  INLINE void operator()(Int i, value_type& update) const {
    update = unite(update, BBox<dim>(get_vec<dim>(coords_, i)));
  }
};

template <Int dim>
BBox<dim> find_bounding_box(Reals coords) {
  CHECK(coords.size() % dim == 0);
  return parallel_reduce(coords.size() / dim, BBoxFunctor<dim>(coords));
}
