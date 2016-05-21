Adj transit(Adj h2m, Adj m2l, Int high_dim, Int low_dim) {
  CHECK(3 >= high_dim);
  Int mid_dim = low_dim + 1;
  CHECK(high_dim > mid_dim);
  CHECK(low_dim == 1 || low_dim == 0);
  LOs hm2m = h2m.ab2b;
  Read<I8> m2hm_codes = h2m.codes;
  LOs ml2l = m2l.ab2b;
  Read<I8> ml_codes = m2l.codes;
  LO nmids_per_high = simplex_degrees[high_dim][mid_dim];
  LO nlows_per_mid = simplex_degrees[mid_dim][low_dim];
  LO nlows_per_high = simplex_degrees[high_dim][low_dim];
  LO nhighs = hm2m.size() / nmids_per_high;
  Write<LO> hl2l(nhighs * nlows_per_high);
  Write<I8> codes;
  if (low_dim == 1)
    codes = Write<I8>(hl2l.size());
  auto f = LAMBDA(LO h) {
    LO hl_begin = h * nlows_per_high;
    LO hm_begin = h * nmids_per_high;
    for (Int hl = 0; hl < nlows_per_high; ++hl) {
      UpTemplate ut = up_templates[high_dim][low_dim][hl][0];
      Int hm = ut.up;
      Int hml = ut.which_down;
      LO m = hm2m[hm_begin + hm];
      auto m2hm_code = m2hm_codes[hm_begin + hm];
      auto hm2m_code = invert_alignment(nlows_per_mid, m2hm_code);
      Int ml = align_index(nlows_per_mid, low_dim, hml, hm2m_code);
      LO ml_begin = m * nlows_per_mid;
      LO l = ml2l[ml_begin + ml];
      //safety check for duplicates.
      //remove after this code is heavily exercised (or don't)
      for (Int hhl2 = 0; hhl2 < hl; ++hhl2)
        CHECK(l != hl2l[hl_begin + hhl2]);
      hl2l[hl_begin + hl] = l;
      if (low_dim == 1) {
        auto tet_tri_code = hm2m_code;
        auto tri_edge_code = ml_codes[ml_begin + ml];
        bool tet_tri_flipped = code_is_flipped(tet_tri_code);
        bool tri_edge_flipped = code_rotation(tri_edge_code) == 1;
        bool canon_flipped = ut.is_flipped;
        bool tet_edge_flipped = tet_tri_flipped ^
                                tri_edge_flipped ^
                                canon_flipped;
        codes[hl_begin + hl] = make_code(false, tet_edge_flipped, 0);
      }
    }
  };
  parallel_for(nhighs, f);
  if (low_dim == 1)
    return Adj(hl2l, codes);
  return Adj(hl2l);
}
