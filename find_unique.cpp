template <Int deg, typename T>
void make_canonical(Read<T> ev2v, Read<T>* canon_out, Read<I8>* codes_out) {
  LO nev = ev2v.size();
  LO ne = nev / deg;
  Write<T> canon(nev);
  Write<I8> codes(ne);
  auto f = LAMBDA(LO e) {
    LO begin = e * deg;
    /* find the smallest vertex */
    Int min_j = 0;
    auto min_v = ev2v[begin];
    for (Int j = 1; j < deg; ++j) {
      LO ev = j + begin;
      auto v = ev2v[ev];
      if (v < min_v) {
        min_j = j;
        min_v = v;
      }
    }
    /* rotate to make it first */
    auto rotation = rotation_to_first<deg>(min_j);
    rotate_adj<deg>(rotation, &ev2v[begin], &canon[begin]);
    auto is_flipped = false;
    if (deg == 3 && canon[begin + 2] < canon[begin + 1]) {
      is_flipped = true;
      flip_adj(&canon[begin]);
    }
    codes[e] = make_code(is_flipped, rotation, 0);
  };
  parallel_for(ne, f);
  *canon_out = canon;
  *codes_out = codes;
}

#define INST_CANON(deg,T) \
template void make_canonical<deg, T>(Read<T> ev2v, \
    Read<T>* canon_out, Read<I8>* codes_out);
INST_CANON(2,LO)
INST_CANON(3,LO)
INST_CANON(2,GO)
INST_CANON(3,GO)
#undef INST_CANON

/* check whether adjacent lists of (deg) vertices
   are the same */
template <Int deg>
DEVICE static bool are_equal(LOs canon, LO e0, LO e1) {
  LO a = e0 * deg;
  LO b = e1 * deg;
  for (LO j = 0; j < deg; ++j)
    if (canon[a + j] != canon[b + j])
      return false;
  return true;
}

template <Int deg>
Read<I8> find_jumps(LOs canon, LOs e_sorted2e) {
  auto ne = e_sorted2e.size();
  Write<I8> jumps(ne, 0);
  auto f = LAMBDA(LO e_sorted) {
    auto e0 = e_sorted2e[e_sorted];
    auto e1 = e_sorted2e[e_sorted + 1];
    if (!are_equal<deg>(canon, e0, e1))
      jumps[e_sorted] = 1;
  };
  parallel_for(ne - 1, f);
  if (jumps.size()) jumps.set(jumps.size() - 1, 1);
  return jumps;
}

template Read<I8> find_jumps<2>(LOs canon, LOs e_sorted2e);
template Read<I8> find_jumps<3>(LOs canon, LOs e_sorted2e);

template <Int deg>
LOs find_unique(LOs uv2v) {
  LOs uv2v_canon;
  Read<I8> u_codes;
  make_canonical<deg>(uv2v, &uv2v_canon, &u_codes);
  auto sorted2u = sort_by_keys<LO,deg>(uv2v_canon);
  auto jumps = find_jumps<deg>(uv2v_canon, sorted2u);
  auto e2sorted = collect_marked(jumps);
  auto e2u = compound_maps(e2sorted, sorted2u);
  return unmap<LO>(e2u, uv2v, deg);
}

LOs find_unique(LOs hv2v, Int high_dim, Int low_dim) {
  auto uv2v = form_uses(hv2v, high_dim, low_dim);
  if (low_dim == 1)
    return find_unique<2>(uv2v);
  if (low_dim == 2)
    return find_unique<3>(uv2v);
  NORETURN(LOs({}));
}
