Adj::Adj() {
}

Adj::Adj(LOs ab2b_):
  ab2b(ab2b_) {
}

Adj::Adj(LOs ab2b_, Read<I8> codes_):
  ab2b(ab2b_),codes(codes_) {
}

Adj::Adj(LOs a2ab_, LOs ab2b_, Read<I8> codes_):
  a2ab(a2ab_),ab2b(ab2b_),codes(codes_) {
}
static void order_by_globals(
    LOs l2lh,
    Write<LO> lh2h,
    Write<I8> codes,
    Read<GO> hg) {
  LO nl = l2lh.size() - 1;
  auto f = LAMBDA(LO l) {
    LO begin = l2lh[l];
    LO end = l2lh[l + 1];
    for (LO j = begin; j < end; ++j) {
      LO k_min = j;
      GO min_g = hg[lh2h[j]];
      for (LO k = j + 1; k < end; ++k) {
        GO g = hg[lh2h[k]];
        if (g < min_g) {
          k_min = k;
          min_g = g;
        }
      }
      swap2(lh2h[j], lh2h[k_min]);
      swap2(codes[j], codes[k_min]);
    }
  };
  parallel_for(nl, f);
}

Adj invert(Adj down, I8 nlows_per_high, LO nlows,
    Read<GO> high_globals) {
  LOs l2lh;
  LOs lh2hl;
  map::invert(down.ab2b, nlows, l2lh, lh2hl, map::BY_ATOMICS);
  LO nlh = lh2hl.size();
  Read<I8> down_codes(down.codes);
  Write<LO> lh2h(nlh);
  Write<I8> codes(nlh);
  if (down_codes.size()) {
    auto f = LAMBDA(LO lh) {
      LO hl = lh2hl[lh];
      LO h = hl / nlows_per_high;
      lh2h[lh] = h;
      I8 which_down = hl % nlows_per_high;
      I8 down_code = down_codes[hl];
      bool is_flipped = code_is_flipped(down_code);
      I8 rotation = code_rotation(down_code);
      codes[lh] = make_code(is_flipped, rotation, which_down);
    };
    parallel_for(nlh, f);
  } else {
    auto f = LAMBDA(LO lh) {
      LO hl = lh2hl[lh];
      LO h = hl / nlows_per_high;
      lh2h[lh] = h;
      I8 which_down = hl % nlows_per_high;
      codes[lh] = make_code(0, 0, which_down);
    };
    parallel_for(nlh, f);
  }
  order_by_globals(l2lh, lh2h, codes, high_globals);
  return Adj(l2lh, lh2h, codes);
}
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
LOs form_uses(LOs hv2v, I8 high_dim, I8 low_dim) {
  I8 nverts_per_high = degrees[high_dim][0];
  I8 nverts_per_low = degrees[low_dim][0];
  I8 nlows_per_high = degrees[high_dim][low_dim];
  LO nhigh = hv2v.size() / nverts_per_high;
  LO nuses = nhigh * nlows_per_high;
  Write<LO> uv2v(nuses * nverts_per_low);
  auto f = LAMBDA(LO h) {
    LO h_begin = h * nverts_per_high;
    for (I8 u = 0; u < nlows_per_high; ++u) {
      LO u_begin = (h * nlows_per_high + u) * nverts_per_low;
      for (I8 uv = 0; uv < nverts_per_low; ++uv) {
        uv2v[u_begin + uv] =
          hv2v[h_begin + simplices[high_dim][low_dim][u][uv]];
      }
    }
  };
  parallel_for(nhigh, f);
  return uv2v;
}
template <Int deg>
void reflect_down(LOs euv2v, LOs ev2v,
    LOs& eu2e, Read<I8>& eu2e_codes_) {
  LO neu = euv2v.size() / deg;
  LOs euv2v_canon;
  Read<I8> eu_codes;
  make_canonical<deg>(euv2v, euv2v_canon, eu_codes);
  LOs ev2v_canon;
  Read<I8> e_codes;
  make_canonical<deg>(ev2v, ev2v_canon, e_codes);
  LOs eu_sorted2eu = sort_by_keys<LO,deg>(euv2v_canon);
  LOs e_sorted2e = sort_by_keys<LO,deg>(ev2v_canon);
  Read<I8> jumps = find_jumps<deg>(euv2v_canon, eu_sorted2eu);
  LOs eu_sorted2e_sorted = excl_scan<LO,I8>(jumps);
  LOs eu2eu_sorted = invert_permutation(eu_sorted2eu);
  LOs eu2e_sorted = compound_maps(eu2eu_sorted, eu_sorted2e_sorted);
  eu2e = compound_maps(eu2e_sorted, e_sorted2e);
  Write<I8> eu2e_codes(neu);
  auto f = LAMBDA(LO eu) {
    LO e = eu2e[eu];
    eu2e_codes[eu] = compound_alignments<deg>(
        e_codes[e], invert_alignment<deg>(eu_codes[eu]));
  };
  parallel_for(neu, f);
  eu2e_codes_ = eu2e_codes;
}

template void reflect_down<2>(LOs euv2v, LOs ev2v,
    LOs& eu2e, Read<I8>& eu2e_codes_);
template void reflect_down<3>(LOs euv2v, LOs ev2v,
    LOs& eu2e, Read<I8>& eu2e_codes_);

Adj reflect_down(LOs hv2v, LOs lv2v, I8 high_dim, I8 low_dim) {
  LOs uv2v = form_uses(hv2v, high_dim, low_dim);
  LOs hl2l;
  Read<I8> codes;
  if (low_dim == 1)
    reflect_down<2>(uv2v, lv2v, hl2l, codes);
  if (low_dim == 2)
    reflect_down<3>(uv2v, lv2v, hl2l, codes);
  return Adj(hl2l, codes);
}
