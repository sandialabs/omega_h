Graph add_edges(Graph g1, Graph g2) {
  auto v2e1 = g1.a2ab;
  auto e2v1 = g1.ab2b;
  auto v2e2 = g2.a2ab;
  auto e2v2 = g2.ab2b;
  auto nv = v2e1.size() - 1;
  auto deg1 = get_degrees(v2e1);
  auto deg2 = get_degrees(v2e2);
  auto deg = add_each(deg1, deg2);
  auto v2e = offset_scan<LO>(deg);
  Write<LO> e2v(v2e.get(v2e.size() - 1));
  auto f = LAMBDA(LO v) {
    auto begin1 = v2e1[v];
    auto end1 = v2e1[v + 1];
    auto begin2 = v2e2[v];
    auto end2 = v2e2[v + 1];
    auto begin = v2e[v];
    auto end = v2e[v + 1];
    auto k = begin;
    for (auto j = begin1; j < end1; ++j)
      e2v[k++] = e2v1[j];
    for (auto j = begin2; j < end2; ++j)
      e2v[k++] = e2v2[j];
    CHECK(k == end);
  };
  parallel_for(nv, f);
  return Graph(v2e, e2v);
}

Graph unmap_graph(LOs a2b, Graph b2c) {
  auto b2bc = b2c.a2ab;
  auto bc2c = b2c.ab2b;
  auto b_degrees = get_degrees(b2bc);
  auto a_degrees = unmap(a2b, b_degrees);
  auto a2ac = offset_scan<LO>(a_degrees);
  auto na = a2b.size();
  Write<LO> ac2c(a2ac.get(a2ac.size() - 1));
  auto f = LAMBDA(LO a) {
    auto b = a2b[a];
    auto bc = b2bc[b];
    for (auto ac = a2ac[a]; ac < a2ac[a + 1]; ++ac) {
      ac2c[ac] = bc2c[bc++];
    }
  };
  parallel_for(na, f);
  return Graph(a2ac, ac2c);
}
