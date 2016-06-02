template <Int sdim, Int edim>
OSH_INLINE Few<Vector<sdim>, edim> simplex_basis(Few<Vector<sdim>, edim + 1> p) {
  Few<Vector<sdim>, edim> b;
  for (Int i = 0; i < edim; ++i)
    b[i] = p[i + 1] - p[0];
  return b;
}

OSH_INLINE Real triangle_area(Few<Vector<2>, 2> b) {
  return cross(b[0], b[1]) / 2.0;
}

OSH_INLINE Real triangle_area(Few<Vector<3>, 2> b) {
  return norm(cross(b[0], b[1])) / 2.0;
}

OSH_INLINE Real tet_volume(Few<Vector<3>, 3> b) {
  return (cross(b[0], b[1]) * b[2]) / 6.0;
}

template <Int dim>
OSH_INLINE Real iso_edge_length(Few<Vector<dim>, 2> p, Real iso) {
  return norm(p[1] - p[0]) / iso;
}

template <Int dim>
OSH_INLINE Real iso_edge_length(Few<LO, 2> v, Reals coords, Reals isos) {
  auto p = gather_vectors<2,dim>(coords, v);
  auto iso = average(gather_scalars<2>(isos, v));
  return iso_edge_length(p, iso);
}

template <Int dim>
OSH_INLINE Real metric_edge_length(Few<Vector<dim>, 2> p, Matrix<dim,dim> metric) {
  return metric_length(metric, p[1] - p[0]);
}

template <Int dim>
OSH_INLINE Real metric_edge_length(Few<LO, 2> v, Reals coords, Reals metrics) {
  auto p = gather_vectors<2,dim>(coords, v);
  auto metric = average_metrics(gather_symms<2,dim>(metrics, v));
  return metric_edge_length(p, metric);
}

template <Int dim>
struct IsoEdgeLengths {
  Reals coords;
  Reals isos;
  IsoEdgeLengths(Mesh const& mesh):
    coords(mesh.coords()),
    isos(mesh.get_array<Real>(VERT, "size"))
  {}
  OSH_INLINE Real measure(Few<LO, 2> v) const {
    return iso_edge_length<dim>(v, coords, isos);
  }
};

template <Int dim>
struct MetricEdgeLengths {
  Reals coords;
  Reals metrics;
  MetricEdgeLengths(Mesh const& mesh):
    coords(mesh.coords()),
    metrics(mesh.get_tag<Real>(VERT, "metric").array())
  {}
  OSH_INLINE Real measure(Few<LO, 2> v) const {
    return metric_edge_length<dim>(v, coords, metrics);
  }
};

Reals measure_edges(Mesh& mesh, LOs ev2v);
Reals measure_edges(Mesh& mesh);
