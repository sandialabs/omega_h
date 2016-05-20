template <Int deg>
void make_canonical(LOs ev2v,
    LOs& canon_, Read<I8>& codes_) {
  LO nev = ev2v.size();
  LO ne = nev / deg;
  Write<LO> canon(nev);
  Write<I8> codes(ne);
  auto f = LAMBDA(LO e) {
    LO begin = e * deg;
    /* find the smallest vertex */
    Int min_j = 0;
    LO min_v = ev2v[begin];
    for (Int j = 1; j < deg; ++j) {
      LO ev = j + begin;
      LO v = ev2v[ev];
      if (v < min_v) {
        min_j = j;
        min_v = v;
      }
    }
    /* rotate to make it first */
    Int rotation = rotation_to_first<deg>(min_j);
    rotate_adj<deg>(rotation, &ev2v[begin], &canon[begin]);
    bool is_flipped = false;
    if (deg == 3 && canon[begin + 2] < canon[begin + 1]) {
      is_flipped = true;
      flip_adj(&canon[begin]);
    }
    codes[e] = make_code(is_flipped, rotation, 0);
  };
  parallel_for(ne, f);
  canon_ = canon;
  codes_ = codes;
}

template void make_canonical<2>(LOs ev2v,
    LOs& canon_, Read<I8>& codes_);
template void make_canonical<3>(LOs ev2v,
    LOs& canon_, Read<I8>& codes_);

/* check whether adjacent lists of (deg) vertices
   are the same */
template <Int deg>
INLINE static bool are_equal(LOs canon, LO e0, LO e1) {
  LO a = e0 * deg;
  LO b = e1 * deg;
  for (LO j = 0; j < deg; ++j)
    if (canon[a + j] != canon[b + j])
      return false;
  return true;
}

template <Int deg>
Read<I8> find_jumps(LOs canon, LOs e_sorted2e) {
  LO ne = e_sorted2e.size();
  Write<I8> jumps(ne, 0);
  auto f = LAMBDA(LO e_sorted) {
    LO e0 = e_sorted2e[e_sorted];
    LO e1 = e_sorted2e[e_sorted + 1];
    if (!are_equal<deg>(canon, e0, e1))
      jumps[e_sorted] = 1;
  };
  parallel_for(ne - 1, f);
  if (jumps.size())
    jumps.set(jumps.size() - 1, 1);
  return jumps;
}

template Read<I8> find_jumps<2>(LOs canon, LOs e_sorted2e);
template Read<I8> find_jumps<3>(LOs canon, LOs e_sorted2e);

template <Int deg>
LOs find_unique(LOs uv2v) {
  LOs uv2v_canon;
  Read<I8> u_codes;
  make_canonical<deg>(uv2v, uv2v_canon, u_codes);
  LOs sorted2u = sort_by_keys<LO,deg>(uv2v_canon);
  Read<I8> jumps = find_jumps<deg>(uv2v_canon, sorted2u);
  LOs e2sorted = collect_marked(jumps);
  LOs e2u = compound_maps(e2sorted, sorted2u);
  return unmap<LO,deg>(e2u, uv2v);
}

LOs find_unique(LOs hv2v, I8 high_dim, I8 low_dim) {
  LOs uv2v = form_uses(hv2v, high_dim, low_dim);
  if (low_dim == 1)
    return find_unique<2>(uv2v);
  if (low_dim == 2)
    return find_unique<3>(uv2v);
  NORETURN(LOs({}));
}
