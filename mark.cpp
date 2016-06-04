Read<I8> mark_exposed_sides(Mesh& mesh) {
  auto ns = mesh.nents(mesh.dim() - 1);
  auto s2sc = mesh.ask_up(mesh.dim() - 1, mesh.dim()).a2ab;
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

template <typename T>
Read<I8> mark_equal(Read<T> a, T val) {
  Write<I8> out(a.size());
  auto f = LAMBDA(LO i) {
    out[i] = (a[i] == val);
  };
  parallel_for(out.size(), f);
  return out;
}

template Read<I8> mark_equal(Read<I8> a, I8 val);
template Read<I8> mark_equal(Read<I32> a, I32 val);
template Read<I8> mark_equal(Read<I64> a, I64 val);
template Read<I8> mark_equal(Read<Real> a, Real val);

Read<I8> mark_by_class_dim(Mesh& mesh, Int ent_dim, Int class_dim) {
  auto e2class_dim = mesh.get_array<I8>(ent_dim, "class_dim");
  return mark_equal(e2class_dim, static_cast<I8>(class_dim));
}
