template <typename T>
Read<T> permute(LOs out2in, Read<T> in) {
  Write<T> out(in.size());
  auto f = LAMBDA(LO i) {
    out[i] = in[out2in[i]];
  };
  return out;
}

template Read<I8  > permute(LOs new2old, Read<I8  > in);
template Read<I32 > permute(LOs new2old, Read<I32 > in);
template Read<I64 > permute(LOs new2old, Read<I64 > in);
template Read<Real> permute(LOs new2old, Read<Real> in);

LOs invert_funnel(LOs ab2a, LO na) {
  LO nab = ab2a.size();
  Write<LO> a2ab(na + 1, -1);
  a2ab.set(0, 0);
  auto f = LAMBDA(LO ab) {
    LO a_end = ab2a[ab];
    LO a_start = ab2a[ab + 1];
    if (a_end != a_start) {
      a2ab[a_end + 1] = ab + 1;
    }
  };
  parallel_for(nab - 1, f);
  if (nab) {
    LO a_end = ab2a.get(nab - 1);
    a2ab.set(a_end + 1, nab);
  }
  fill_right(a2ab);
  return a2ab;
}

void invert_map_by_sort(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a) {
  LOs ab2b = a2b;
  LOs ba2ab = sort_by_keys(ab2b);
  LOs ba2b = permute(ba2ab, ab2b);
  b2ba = invert_funnel(ba2b, nb);
  ba2a = ba2ab;
}
