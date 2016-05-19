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
