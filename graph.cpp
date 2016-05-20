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
