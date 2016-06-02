template <typename EdgeLengths>
Reals measure_edges_tmpl(Mesh& mesh, LOs ev2v) {
  EdgeLengths measurer(mesh);
  auto nedges = ev2v.size() / 2;
  Write<Real> lengths(nedges);
  auto f = OSH_LAMBDA(LO e) {
    auto v = gather_verts<2>(ev2v, e);
    lengths[e] = measurer.measure(v);
  };
  parallel_for(nedges, f);
  return lengths;
}

Reals measure_edges(Mesh& mesh, LOs ev2v) {
  if (mesh.dim() == 3) {
    if (mesh.has_tag(VERT, "size")) {
      return measure_edges_tmpl<IsoEdgeLengths<3>>(mesh, ev2v);
    }
    if (mesh.has_tag(VERT, "metric")) {
      return measure_edges_tmpl<MetricEdgeLengths<3>>(mesh, ev2v);
    }
  } else {
    CHECK(mesh.dim() == 2);
    if (mesh.has_tag(VERT, "size")) {
      return measure_edges_tmpl<IsoEdgeLengths<2>>(mesh, ev2v);
    }
    if (mesh.has_tag(VERT, "metric")) {
      return measure_edges_tmpl<MetricEdgeLengths<2>>(mesh, ev2v);
    }
  }
  fail("measure_edges(): no size field exists!\n");
}

Reals measure_edges(Mesh& mesh) {
  return measure_edges(mesh, mesh.ask_verts_of(EDGE));
}
