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
      Int which_down = hl % nlows_per_high;
      I8 down_code = down_codes[hl];
      bool is_flipped = code_is_flipped(down_code);
      Int rotation = code_rotation(down_code);
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
