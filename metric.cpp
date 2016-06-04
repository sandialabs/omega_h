template <Int sdim, Int edim>
Reals average_metric_tmpl(Mesh& mesh, LOs entities) {
  auto nents = entities.size();
  Write<Real> out(nents * symm_dofs(sdim));
  auto ev2v = mesh.ask_verts_of(edim);
  auto v_metrics = mesh.get_array<Real>(VERT, "metric");
  auto f = LAMBDA(LO i) {
    auto e = entities[i];
    auto v = gather_verts<edim + 1>(ev2v, e);
    auto ms = gather_symms<edim + 1, sdim>(v_metrics, v);
    auto m = average_metrics(ms);
    set_symm(out, i, m);
  };
  parallel_for(nents, f);
  return out;
}

Reals average_metric(Mesh& mesh, Int ent_dim, LOs entities) {
  if (mesh.dim() == 3) {
    if (ent_dim == 3) {
      return average_metric_tmpl<3,3>(mesh, entities);
    } else if (ent_dim == 2) {
      return average_metric_tmpl<3,2>(mesh, entities);
    } else {
      CHECK(ent_dim == 1);
      return average_metric_tmpl<3,1>(mesh, entities);
    }
  } else {
    CHECK(mesh.dim() == 2);
    if (ent_dim == 2) {
      return average_metric_tmpl<2,2>(mesh, entities);
    } else {
      CHECK(ent_dim == 1);
      return average_metric_tmpl<2,1>(mesh, entities);
    }
  }
}
