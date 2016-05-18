/* for each entity (or entity use), sort its vertex list.
   express the sorting transformation as
   an alignment code, output those too */
template <Int deg>
static void make_canonical(LOs ev2v,
    LOs& canon_, Read<I8>& codes_) {
  LO nev = ev2v.size();
  LO ne = nev / deg;
  Write<LO> canon(nev);
  Write<I8> codes(ne);
  auto f = LAMBDA(LO e) {
    LO begin = e * deg;
    /* find the smallest vertex */
    I8 min_j = 0;
    LO min_v = ev2v[begin];
    for (I8 j = 1; j < deg; ++j) {
      LO ev = j + begin;
      LO v = ev2v[ev];
      if (v < min_v) {
        min_j = j;
        min_v = v;
      }
    }
    /* rotate to make it first */
    I8 rotation = rotation_to_first<deg>(min_j);
    rotate_adj<deg>(rotation, &ev2v[begin], &canon[begin]);
    bool is_flipped = false;
    if (deg == 3 && canon[2] < canon[1]) {
      is_flipped = true;
      flip_adj(&canon[begin]);
    }
    codes[e] = make_code(is_flipped, rotation, 0);
  };
  parallel_for(ne, f);
  canon_ = canon;
  codes_ = codes;
}

/* check whether adjacent lists of (deg) vertices
   are the same */
template <Int deg>
static bool is_jump(LOs canon, LO h) {
  LO a = h * deg;
  LO b = (h + 1) * deg;
  for (LO j = 0; j < deg; ++j)
    if (canon[a + j] != canon[b + j])
      return true;
  return false;
}

template <Int deg>
static Read<I8> find_jumps(LOs canon) {
  LO nh = canon.size() / deg;
  Write<I8> jumps(nh, 0);
  auto f = LAMBDA(LO h) {
    jumps[h] = is_jump<deg>(canon, h);
  };
  parallel_for(nh - 1, f);
  return jumps;
}

/* given entity uses and unique entities,
   both defined by vertex lists, match
   uses to unique entities including their
   respective alignment codes */
template <Int deg>
void reflect_down(LOs euv2v, LOs ev2v,
    LOs& eu2e, Read<I8>& eu2e_codes_) {
  LO neu = euv2v.size() / deg;
  LOs euv2v_canon;
  Read<I8> eu_codes;
  make_canonical<deg>(euv2v, euv2v_canon, eu_codes);
  std::cerr << "euv2v_canon\n" << euv2v_canon << '\n';
  LOs ev2v_canon;
  Read<I8> e_codes;
  make_canonical<deg>(ev2v, ev2v_canon, e_codes);
  LOs eu_sorted2eu = sort_by_keys<LO,deg>(euv2v_canon);
  std::cerr << "eu_sorted2eu\n" << eu_sorted2eu << '\n';
  LOs e_sorted2e = sort_by_keys<LO,deg>(ev2v_canon);
  std::cerr << "e_sorted2e\n" << e_sorted2e << '\n';
  Read<I8> jumps = find_jumps<deg>(euv2v_canon);
  std::cerr << "jumps\n" << jumps << '\n';
  LOs eu_sorted2e_sorted = excl_scan<LO,I8>(jumps);
  std::cerr << "eu_sorted2e_sorted\n" << eu_sorted2e_sorted << '\n';
  LOs eu2eu_sorted = invert_permutation(eu_sorted2eu);
  std::cerr << "eu2eu_sorted\n" << eu2eu_sorted << '\n';
  LOs eu2e_sorted = compound_maps(eu2eu_sorted, eu_sorted2e_sorted);
  std::cerr << "eu2e_sorted\n" << eu2e_sorted << '\n';
  eu2e = compound_maps(eu2e_sorted, e_sorted2e);
  Write<I8> eu2e_codes(neu);
  auto f = LAMBDA(LO eu) {
    LO e = eu2e[eu];
    eu2e_codes[eu] = compound_alignments<deg>(
        eu_codes[eu], invert_alignment<deg>(e_codes[e]));
  };
  parallel_for(neu, f);
  eu2e_codes_ = eu2e_codes;
}

template void reflect_down<2>(LOs euv2v, LOs ev2v,
    LOs& eu2e, Read<I8>& eu2e_codes_);
template void reflect_down<3>(LOs euv2v, LOs ev2v,
    LOs& eu2e, Read<I8>& eu2e_codes_);
