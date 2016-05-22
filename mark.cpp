Read<I8> mark_exposed_sides(Mesh& mesh) {
  auto ns = mesh.nents(mesh.dim() - 1);
  auto s2sc = mesh.ask_adj(mesh.dim() - 1, mesh.dim()).a2ab;
  Write<I8> exposed(ns);
  auto f = LAMBDA(LO s) {
    exposed[s] = ((s2sc[s + 1] - s2sc[s]) < 2);
  };
  parallel_for(ns, f);
  return exposed;
}

Read<I8> mark_down(Mesh& mesh, Int high_dim, Int low_dim,
    Read<I8> high_marked) {
  auto l2h = mesh.ask_up(low_dim, high_dim);
  auto l2lh = l2h.a2ab;
  auto lh2h = l2h.ab2b;
  auto nl = mesh.nents(low_dim);
  Write<I8> out(nl, 0);
  auto f = LAMBDA(LO l) {
    for (LO lh = l2lh[l]; lh < l2lh[l + 1]; ++lh)
      if (high_marked[lh2h[lh]])
        out[l] = 1;
  };
  parallel_for(nl, f);
  return out;
}

Read<I8> mark_by_class_dim(Mesh& mesh, Int ent_dim, Int class_dim) {
  auto ne = mesh.nents(ent_dim);
  auto e2class_dim = mesh.get_tag<I8>(ent_dim, "class_dim").array();
  Write<I8> out(ne);
  auto f = LAMBDA(LO e) {
    out[e] = (e2class_dim[e] == class_dim);
  };
  parallel_for(ne, f);
  return out;
}
