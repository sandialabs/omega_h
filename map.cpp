template <typename T>
Read<T> permute(LOs out2in, Read<T> in) {
  Write<T> out(in.size());
  auto f = LAMBDA(LO i) {
    out[i] = in[out2in[i]];
  };
  parallel_for(out.size(), f);
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

namespace map {

void invert_by_sorting(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a) {
  LOs ab2b = a2b;
  LOs ba2ab = sort_by_keys(ab2b);
  LOs ba2b = permute(ba2ab, ab2b);
  b2ba = invert_funnel(ba2b, nb);
  ba2a = ba2ab;
}

void invert_by_atomics(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a) {
  LO na = a2b.size();
  Write<LO> degrees(nb, 0);
  auto count = LAMBDA(LO a) {
    atomic_increment(&degrees[a2b[a]]);
  };
  parallel_for(na, count);
  b2ba = offset_scan<LO>(Read<LO>(degrees));
  LO nba = b2ba.get(nb);
  Write<LO> write_ba2a(nba);
  degrees = Write<LO>(nb, 0);
  auto fill = LAMBDA(LO a) {
    LO b = a2b[a];
    LO first = b2ba[b];
    LO j = atomic_fetch_add<LO>(&degrees[a2b[a]], 1);
    write_ba2a[first + j] = a;
  };
  parallel_for(na, fill);
  ba2a = write_ba2a;
}

void invert(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a, InvertMethod method) {
  if (method == BY_SORTING)
    invert_by_sorting(a2b, nb, b2ba, ba2a);
  if (method == BY_ATOMICS)
    invert_by_atomics(a2b, nb, b2ba, ba2a);
}

}
