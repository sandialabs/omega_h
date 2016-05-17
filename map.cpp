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
    LO a_end = ab2a[nab - 1];
    a2ab.set(a_end + 1, nab);
  }
  fill_right(a2ab);
  return a2ab;
}

void invert_map_by_sorting(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a) {
  LOs ab2b = a2b;
  LOs ba2ab = sort_by_keys(ab2b);
  LOs ba2b = permute(ba2ab, ab2b);
  b2ba = invert_funnel(ba2b, nb);
  ba2a = ba2ab;
}

void invert_map_by_atomics(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a) {
  LO na = a2b.size();
  Write<I8> degrees(nb, 0);
  auto count = LAMBDA(LO a) {
    atomic_increment(&degrees[a2b[a]]);
  };
  parallel_for(na, count);
  b2ba = offset_scan<LO>(Read<I8>(degrees));
  LO nba = b2ba.get(nb);
  Write<LO> write_ba2a(nba);
  degrees = Write<I8>(nb, 0);
  auto fill = LAMBDA(LO a) {
    LO b = a2b[a];
    LO first = b2ba[b];
    LO j = atomic_fetch_add<I8>(&degrees[a2b[a]], 1);
    write_ba2a[first + j] = a;
  };
  parallel_for(na, fill);
  ba2a = write_ba2a;
}

LOs order_by_globals(
    LOs a2ab, LOs ab2b, Read<GO> b_global) {
  LO na = a2ab.size();
  LO nab = ab2b.size();
  Write<LO> write_ab2b(nab);
  auto f = LAMBDA(LO a) {
    LO begin = a2ab[a];
    LO end = a2ab[a + 1];
    for (LO j = begin; j < end; ++j)
      write_ab2b[j] = ab2b[j];
    for (LO j = begin; j < end; ++j) {
      LO k_min = j;
      GO min_g = b_global[write_ab2b[j]];
      for (LO k = j + 1; k < end; ++k) {
        LO b = write_ab2b[k];
        GO g = b_global[b];
        if (g < min_g) {
          k_min = k;
          min_g = g;
        }
      }
      swap2(write_ab2b[j], write_ab2b[k_min]);
    }
  };
  parallel_for(na, f);
  return write_ab2b;
}
