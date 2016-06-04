template <typename ElementQualities, Int dim>
Reals measure_qualities_tmpl(Mesh const& mesh, LOs ev2v) {
  ElementQualities measurer(mesh);
  auto neev = dim + 1;
  auto nelems = ev2v.size() / neev;
  Write<Real> qualities(nelems);
  auto f = LAMBDA(LO e) {
    auto v = gather_verts<dim + 1>(ev2v, e);
    qualities[e] = measurer.measure(v);
  };
  parallel_for(nelems, f);
  return qualities;
}

Reals measure_qualities(Mesh const& mesh, LOs ev2v) {
  if (mesh.dim() == 3) {
    if (mesh.has_tag(VERT, "metric")) {
      return measure_qualities_tmpl<MetricElementQualities,3>(mesh, ev2v);
    } else {
      return measure_qualities_tmpl<RealElementQualities,3>(mesh, ev2v);
    }
  } else {
    CHECK(mesh.dim() == 2);
    if (mesh.has_tag(VERT, "metric")) {
      return measure_qualities_tmpl<MetricElementQualities,2>(mesh, ev2v);
    } else {
      return measure_qualities_tmpl<RealElementQualities,2>(mesh, ev2v);
    }
  }
}

Reals measure_qualities(Mesh& mesh) {
  return measure_qualities(mesh, mesh.ask_verts_of(mesh.dim()));
}
