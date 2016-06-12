template <typename EdgeLengths>
Reals measure_edges_tmpl(Mesh* mesh, LOs a2e) {
  EdgeLengths measurer(mesh);
  auto ev2v = mesh.ask_verts_of(EDGE);
  auto na = a2e.size();
  Write<Real> lengths(na);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<2>(ev2v, e);
    lengths[a] = measurer.measure(v);
  };
  parallel_for(na, f);
  return lengths;
}

Reals measure_edges(Mesh* mesh, LOs a2e) {
  if (mesh.dim() == 3) {
    if (mesh.has_tag(VERT, "size")) {
      return measure_edges_tmpl<IsoEdgeLengths<3>>(mesh, a2e);
    }
    if (mesh.has_tag(VERT, "metric")) {
      return measure_edges_tmpl<MetricEdgeLengths<3>>(mesh, a2e);
    }
  } else {
    CHECK(mesh.dim() == 2);
    if (mesh.has_tag(VERT, "size")) {
      return measure_edges_tmpl<IsoEdgeLengths<2>>(mesh, a2e);
    }
    if (mesh.has_tag(VERT, "metric")) {
      return measure_edges_tmpl<MetricEdgeLengths<2>>(mesh, a2e);
    }
  }
  fail("measure_edges(): no size field exists!\n");
}

Reals measure_edges(Mesh* mesh) {
  return measure_edges(mesh, LOs(mesh.nedges(), 0, 1));
}
