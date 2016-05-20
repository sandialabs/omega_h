Adj transit(Adj h2m, Adj m2l, I8 high_dim, I8 low_dim) {
  std::cerr << "transit high_dim " << Int(high_dim) << " low_dim " << Int(low_dim) << '\n';
  CHECK(3 >= high_dim);
  I8 mid_dim = low_dim + 1;
  CHECK(high_dim > mid_dim);
  CHECK(low_dim == 1 || low_dim == 0);
  LOs hm2m = h2m.ab2b;
  Read<I8> m2hm_codes = h2m.codes;
  LOs ml2l = m2l.ab2b;
  Read<I8> ml_codes = m2l.codes;
  LO nmids_per_high = degrees[high_dim][mid_dim];
  LO nlows_per_mid = degrees[mid_dim][low_dim];
  LO nlows_per_high = degrees[high_dim][low_dim];
  LO nhighs = hm2m.size() / nmids_per_high;
  Write<LO> hl2l(nhighs * nlows_per_high);
  Write<I8> codes;
  if (low_dim == 1)
    codes = Write<I8>(hl2l.size());
  std::cerr << "hm2m\n" << hm2m << '\n';
  std::cerr << "ml2l\n" << ml2l << '\n';
  auto f = LAMBDA(LO h) {
    LO hl_begin = h * nlows_per_high;
    LO hm_begin = h * nmids_per_high;
    for (LO hl = 0; hl < nlows_per_high; ++hl) {
      UpTemplate ut = up_templates[high_dim][low_dim][hl][0];
      std::cerr << singular_names[low_dim] << " " << hl
        << " of a " << singular_names[high_dim] << "...\n";
      Int hm = ut.up;
      Int hml = ut.which_down;
      std::cerr << "is the " << hml << "'th " << singular_names[low_dim]
        << " of the the ";
      std::cerr << singular_names[high_dim] << "'s " << hm
        << "'th " << singular_names[mid_dim] << '\n';
      LO m = hm2m[hm_begin + hm];
      std::cerr << "for " << singular_names[high_dim] << " " << h
        << ", that is " << singular_names[mid_dim] << " " << m << '\n';
      I8 m2hm_code = m2hm_codes[hm_begin + hm];
      std::cerr << "which it is using with rotation " << code_rotation(m2hm_code);
      if (code_is_flipped(m2hm_code))
        std::cerr << " and a flip";
      std::cerr << '\n';
      I8 hm2m_code = invert_alignment(nlows_per_mid, m2hm_code);
      std::cerr << "getting from the " << singular_names[mid_dim]
        << " use to " << singular_names[mid_dim] << " " << m
        << " would take rotation " << code_rotation(hm2m_code);
      if (code_is_flipped(hm2m_code))
        std::cerr << " and a flip";
      std::cerr << '\n';
      Int ml = align_index(nlows_per_mid, hml, hm2m_code);
      std::cerr << "making the " << hml << "'th " << singular_names[low_dim]
        << " actually the " << ml << "'th\n";
      LO ml_begin = m * nlows_per_mid;
      LO l = ml2l[ml_begin + ml];
      std::cerr << "which is " << singular_names[low_dim] << " " << l << '\n';
      hl2l[hl_begin + hl] = l;
      if (low_dim == 1) {
        I8 tet_tri_code = hm2m_code;
        I8 tri_edge_code = ml_codes[ml_begin + ml];
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
  std::cerr << "hl2l\n" << LOs(hl2l) << '\n';
  if (low_dim == 1)
    return Adj(hl2l, codes);
  return Adj(hl2l);
}
