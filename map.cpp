template <typename T>
void map_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width) {
  auto na = a2b.size();
  CHECK(a_data.size() == na * width);
  auto f = LAMBDA(LO a) {
    auto b = a2b[a];
    for (Int j = 0; j < width; ++j) {
      b_data[b * width + j] = a_data[a * width + j];
    }
  };
  parallel_for(na, f);
}

template <typename T>
Read<T> unmap(LOs a2b, Read<T> b_data, Int width) {
  auto na = a2b.size();
  Write<T> a_data(na * width);
  auto f = LAMBDA(LO a) {
    auto b = a2b[a];
    for (Int j = 0; j < width; ++j) {
      a_data[a * width + j] = b_data[b * width + j];
    }
  };
  parallel_for(na, f);
  return a_data;
}

template <typename T>
Read<T> expand(Read<T> a_data, LOs a2b, Int width) {
  auto na = a2b.size() - 1;
  auto nb = a2b.last();
  CHECK(a_data.size() == na * width);
  Write<T> b_data(nb * width);
  auto f = LAMBDA(LO a) {
    for (auto b = a2b[a]; b < a2b[a + 1]; ++b) {
      for (Int j = 0; j < width; ++j) {
        b_data[b * width + j] = a_data[a * width + j];
      }
    }
  };
  parallel_for(na, f);
  return b_data;
}

template <typename T>
Read<T> permute(Read<T> a_data, LOs a2b, Int width) {
  auto na = a2b.size();
  auto nb = na;
  Write<T> b_data(nb * width);
  auto f = LAMBDA(LO a) {
    auto b = a2b[a];
    for (Int j = 0; j < width; ++j) {
      b_data[b * width + j] = a_data[a * width + j];
    }
  };
  parallel_for(na, f);
  return b_data;
}

LOs multiply_fans(LOs a2b, LOs a2c) {
  auto b_degrees = get_degrees(a2b);
  auto c_degrees = get_degrees(a2c);
  auto degrees = multiply_each(b_degrees, c_degrees);
  auto a2bc = offset_scan(degrees);
  return a2bc;
}

#define INST_T(T) \
template void map_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width); \
template Read<T> unmap(LOs a2b, Read<T> b_data, Int width); \
template Read<T> expand(Read<T> a_data, LOs a2b, Int width); \
template Read<T> permute(Read<T> a_data, LOs a2b, Int width);
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T

LOs compound_maps(LOs a2b, LOs b2c) {
  LO na = a2b.size();
  Write<LO> a2c(a2b.size());
  auto f = LAMBDA(LO a) {
    LO b = a2b[a];
    LO c = b2c[b];
    a2c[a] = c;
  };
  parallel_for(na, f);
  return a2c;
}

LOs invert_permutation(LOs a2b) {
  return invert_injective_map(a2b, a2b.size());
}

Read<I8> invert_marks(Read<I8> marks) {
  Write<I8> out(marks.size());
  auto f = LAMBDA(LO i) {
    out[i] = !marks[i];
  };
  parallel_for(out.size(), f);
  return out;
}

LOs collect_marked(Read<I8> marks) {
  auto ntotal = marks.size();
  auto offsets = offset_scan(marks);
  auto nmarked = offsets.last();
  Write<LO> marked(nmarked);
  auto f = LAMBDA(LO i) {
    if (marks[i])
      marked[offsets[i]] = i;
  };
  parallel_for(ntotal, f);
  return marked;
}

LOs invert_injective_map(LOs a2b, LO nb) {
  Write<LO> b2a(nb, -1);
  auto f = LAMBDA(LO a) {
    b2a[a2b[a]] = a;
  };
  parallel_for(a2b.size(), f);
  return b2a;
}

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
  LOs ba2b = unmap(ba2ab, ab2b, 1);
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
  b2ba = offset_scan(Read<LO>(degrees));
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

LOs get_degrees(LOs offsets) {
  Write<LO> degrees(offsets.size() - 1);
  auto f = LAMBDA(LO i) {
    degrees[i] = offsets[i + 1] - offsets[i];
  };
  parallel_for(degrees.size(), f);
  return degrees;
}

LOs invert_fan(LOs a2b) {
  auto na = a2b.size() - 1;
  auto nb = a2b.last();
  Write<LO> b2a(nb, -1);
  auto f = LAMBDA(LO a)
  {
    if (a2b[a] != a2b[a + 1])
      b2a[a2b[a]] = a;
  };
  parallel_for(na, f);
  fill_right(b2a);
  return b2a;
}
