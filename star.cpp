Adj verts_across_edges(Adj e2v, Adj v2e) {
  auto ev2v = e2v.ab2b;
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ve2e_codes = v2e.codes;
  auto v2vv = v2ve;
  Write<LO> vv2v(ve2e.size());
  auto f = LAMBDA(LO ve) {
    auto vv = ve;
    auto e = ve2e[ve];
    auto ve2e_code = ve2e_codes[ve];
    auto eev = code_which_down(ve2e_code);
    auto v = ev2v[e * 2 + (1 - eev)];
    vv2v[vv] = v;
  };
  parallel_for(vv2v.size(), f);
  return Adj(v2vv, vv2v);
}

Adj edges_across_tris(Adj f2e, Adj e2f) {
  auto fe2e = f2e.ab2b;
  auto e2ef = e2f.a2ab;
  auto ef2f = e2f.ab2b;
  auto ef2f_codes = e2f.codes;
  auto ne = e2ef.size() - 1;
  auto e2ef_degrees = get_degrees(e2ef);
  auto e2ee_degrees = multiply_each_by(2, e2ef_degrees);
  auto e2ee = offset_scan<LO>(e2ee_degrees);
  auto nee = e2ee.get(e2ee.size() - 1);
  Write<LO> ee2e(nee);
  auto lambda = LAMBDA(LO e) {
    auto ef_begin = e2ef[e];
    auto ef_end = e2ef[e + 1];
    auto neef = ef_end - ef_begin;
    auto ee_begin = e2ee[e];
    for (Int eef = 0; eef < neef; ++eef) {
      auto ef = ef_begin + eef;
      auto f = ef2f[ef];
      auto ef2f_code = ef2f_codes[ef];
      auto ffe = code_which_down(ef2f_code);
      auto e1 = fe2e[f * 3 + ((ffe + 1) % 3)];
      auto e2 = fe2e[f * 3 + ((ffe + 2) % 3)];
      ee2e[ee_begin + ef * 2 + 0] = e1;
      ee2e[ee_begin + ef * 2 + 1] = e2;
    }
  };
  parallel_for(ne, lambda);
  return Adj(e2ee, ee2e);
}
