template <typename T>
Read<T> permute(LOs out2in, Read<T> in) {
  Write<T> out(in.size());
  auto f = LAMBDA(LO i) {
    out[i] = in[out2in[i]];
  };
  return out;
}

LOs invert_funnel(LOs ab2a, LO na) {
  LO nab = ab2a.size();
  Write<LO> a2ab(na + 1, -1);
  a2ab.set(0, 0);
  a2ab.set(na, nab);
  auto f = LAMBDA(LO ab) {
    LO a_end = ab2a[ab];
    LO a_start = ab2a[ab + 1];
    if (a_end != a_start) {
      a2ab[a_end + 1] = ab + 1;
    }
  };
  parallel_for(nab - 1, f);
  fill_right(a2ab);
  return a2ab;
}

template Read<I8> permute(LOs new2old, Read<I8> in);
template Read<I32> permute(LOs new2old, Read<I32> in);
template Read<I64> permute(LOs new2old, Read<I64> in);
template Read<Real> permute(LOs new2old, Read<Real> in);
