template <Int dim>
struct BBox {
  OSH_INLINE BBox() {}
  OSH_INLINE BBox(Vector<dim> x):min(x),max(x) {}
  OSH_INLINE BBox(Vector<dim> min_, Vector<dim> max_):min(min_),max(max_) {}
  Vector<dim> min;
  Vector<dim> max;
  /* playing the volatile game again (see int128.hpp) */
  OSH_INLINE void operator=(BBox<dim> const& rhs) volatile {
    min = rhs.min;
    max = rhs.max;
  }
  OSH_INLINE BBox(BBox<dim> const& rhs):
    min(rhs.min),max(rhs.max) {
  }
  OSH_INLINE BBox(const volatile BBox<dim>& rhs):
    min(rhs.min),max(rhs.max) {
  }
};

template <Int dim>
OSH_INLINE BBox<dim> unite(BBox<dim> a, BBox<dim> b) {
  BBox<dim> c;
  for (Int i = 0; i < dim; ++i) {
    c.min[i] = min2(a.min[i], b.min[i]);
    c.max[i] = max2(a.max[i], b.max[i]);
  }
  return c;
}

template <Int dim>
OSH_INLINE bool are_close(BBox<dim> a, BBox<dim> b) {
  return are_close(a.min, b.min) && are_close(a.max, b.max);
}

template <Int dim>
BBox<dim> find_bounding_box(Reals coords);
