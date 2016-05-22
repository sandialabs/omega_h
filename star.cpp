Graph verts_across_edges(Adj e2v, Adj v2e) {
  auto ev2v = e2v.ab2b;
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto v2ve_codes = v2e.codes;
  auto v2vv = v2ve;
  Write<LO> vv2v(ve2e.size());
  auto f = LAMBDA(LO ve) {
    auto vv = ve;
    auto e = ve2e[ve];
    auto v2ve_code = v2ve_codes[ve];
    auto eev = code_which_down(v2ve_code);
    auto v = ev2v[e * 2 + (1 - eev)];
    vv2v[vv] = v;
  };
  parallel_for(vv2v.size(), f);
  return Adj(v2vv, vv2v);
}

Graph edges_across_tris(Adj f2e, Adj e2f) {
  auto fe2e = f2e.ab2b;
  auto e2ef = e2f.a2ab;
  auto ef2f = e2f.ab2b;
  auto e2ef_codes = e2f.codes;
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
      auto e2ef_code = e2ef_codes[ef];
      auto ffe = code_which_down(e2ef_code);
      auto e1 = fe2e[f * 3 + ((ffe + 1) % 3)];
      auto e2 = fe2e[f * 3 + ((ffe + 2) % 3)];
      ee2e[ee_begin + eef * 2 + 0] = e1;
      ee2e[ee_begin + eef * 2 + 1] = e2;
    }
  };
  parallel_for(ne, lambda);
  return Adj(e2ee, ee2e);
}

Graph edges_across_tets(Adj r2e, Adj e2r) {
  auto re2e = r2e.ab2b;
  auto e2er = e2r.a2ab;
  auto er2r = e2r.ab2b;
  auto e2er_codes = e2r.codes;
  auto ne = e2er.size() - 1;
  auto e2ee = e2er;
  Write<LO> ee2e(er2r.size());
  auto f = LAMBDA(LO e) {
    auto er_begin = e2er[e];
    auto er_end = e2er[e + 1];
    for (auto er = er_begin; er < er_end; ++ er) {
      auto r = er2r[er];
      auto e2er_code = e2er_codes[er];
      auto rre = code_which_down(e2er_code);
      auto rre_opp = tet_edges_opp_edges[rre];
      auto re_begin = r * 6;
      auto e_opp = re2e[re_begin + rre_opp];
      auto ee = er;
      ee2e[ee] = e_opp;
    }
  };
  parallel_for(ne, f);
  return Adj(e2ee, ee2e);
}

Graph elements_across_sides(Int dim,
    Adj elems2sides, Adj sides2elems,
    Read<I8> side_is_exposed) {
  auto elem_side2side = elems2sides.ab2b;
  auto side2side_elems = sides2elems.a2ab;
  auto side_elem2elem = sides2elems.ab2b;
  Int nsides_per_elem = dim + 1;
  auto nelems = elem_side2side.size() / nsides_per_elem;
  Write<LO> degrees(nelems);
  auto count = LAMBDA(LO elem) {
    auto begin = elem * nsides_per_elem;
    auto end = begin + nsides_per_elem;
    Int n = 0;
    for (auto elem_side = begin; elem_side < end; ++elem_side) {
      auto side = elem_side2side[elem_side];
      n += !side_is_exposed[side];
    }
    degrees[elem] = n;
  };
  parallel_for(nelems, count);
  auto elem2elem_elems = offset_scan<LO,LO>(degrees);
  auto nelem_elems = elem2elem_elems.get(elem2elem_elems.size() - 1);
  Write<LO> elem_elem2elem(nelem_elems);
  auto fill = LAMBDA(LO elem) {
    auto begin = elem * nsides_per_elem;
    auto end = begin + nsides_per_elem;
    LO elem_elem = elem2elem_elems[elem];
    for (auto elem_side = begin; elem_side < end; ++elem_side) {
      auto side = elem_side2side[elem_side];
      if (!side_is_exposed[side]) {
        auto side_elem = side2side_elems[side];
        if (side_elem2elem[side_elem] == elem)
          ++side_elem;
        elem_elem2elem[elem_elem] = side_elem2elem[side_elem];
        ++elem_elem;
      }
    }
  };
  parallel_for(nelems, fill);
  return Graph(elem2elem_elems, elem_elem2elem);
}
