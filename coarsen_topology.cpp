LOs coarsen_topology(Mesh& mesh,
    LOs keys2verts_onto,
    Int dom_dim,
    Adj keys2doms) {
  auto nccv = simplex_degrees[dom_dim][VERT];
  auto cv2v = mesh.ask_verts_of(dom_dim);
  auto k2kc = keys2doms.a2ab;
  auto kc2c = keys2doms.ab2b;
  auto kc_codes = keys2doms.codes;
  auto ndoms = k2kc.last();
  auto nprods = ndoms;
  auto prod_verts2verts = Write<LO>(nprods * nccv);
  auto nkeys = keys2verts_onto.size();
  auto f = LAMBDA(LO key) {
    auto v_onto = keys2verts_onto[key];
    for (auto kc = k2kc[key]; kc < k2kc[key + 1]; ++kc) {
      auto prod = kc;
      auto c = kc2c[kc];
      auto kc_code = kc_codes[kc];
      auto ccv_col = code_which_down(kc_code);
      auto ppv2v = &prod_verts2verts[prod * nccv];
      for (Int ccv = 0; ccv < nccv; ++ccv) {
        ppv2v[ccv] = cv2v[c * nccv + ccv];
      }
      ppv2v[ccv_col] = v_onto;
    }
  };
  parallel_for(nkeys, f);
  return prod_verts2verts;
}
