template <Int deg>
LOs find_unique(LOs uv2v, I8 nuses_per_high) {
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
  LO nuses_per_high = degrees[high_dim][low_dim];
  if (low_dim == 1)
    return find_unique<2>(uv2v, nuses_per_high);
  if (low_dim == 2)
    return find_unique<3>(uv2v, nuses_per_high);
  NORETURN(LOs({}));
}
