#include "Omega_h_map.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_atomics.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_sort.hpp"

namespace Omega_h {

template <typename T>
void add_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width) {
  auto na = a2b.size();
  OMEGA_H_CHECK(a_data.size() == na * width);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto b = a2b[a];
    for (Int j = 0; j < width; ++j) {
      b_data[b * width + j] += a_data[a * width + j];
    }
  };
  parallel_for(na, f, "add_into");
}

template <typename T>
void map_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width) {
  auto na = a2b.size();
  OMEGA_H_CHECK(a_data.size() == na * width);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto b = a2b[a];
    for (Int j = 0; j < width; ++j) {
      b_data[b * width + j] = a_data[a * width + j];
    }
  };
  parallel_for(na, f, "map_into");
}

template <typename T>
void map_into_range(
    Read<T> a_data, LO begin, LO end, Write<T> b_data, Int width) {
  auto na = end - begin;
  OMEGA_H_CHECK(a_data.size() == na * width);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto b = begin + a;
    for (Int j = 0; j < width; ++j) {
      b_data[b * width + j] = a_data[a * width + j];
    }
  };
  parallel_for(na, f, "map_into_range");
}

template <typename T>
Read<T> map_onto(Read<T> a_data, LOs a2b, LO nb, T init_val, Int width) {
  auto out = Write<T>(nb * width, init_val);
  map_into(a_data, a2b, out, width);
  return out;
}

template <typename T>
Write<T> unmap(LOs a2b, Read<T> b_data, Int width) {
  auto na = a2b.size();
  Write<T> a_data(na * width);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto b = a2b[a];
    for (Int j = 0; j < width; ++j) {
      a_data[a * width + j] = b_data[b * width + j];
    }
  };
  parallel_for(na, f, "unmap");
  return a_data;
}

template <typename T>
Read<T> unmap_range(LO begin, LO end, Read<T> b_data, Int width) {
  OMEGA_H_CHECK(begin <= end);
  auto na = end - begin;
  Write<T> a_data(na * width);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto b = a + begin;
    for (Int j = 0; j < width; ++j) {
      a_data[a * width + j] = b_data[b * width + j];
    }
  };
  parallel_for(na, f, "unmap_range");
  return a_data;
}

template <typename T>
void expand_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width) {
  auto na = a2b.size() - 1;
  OMEGA_H_CHECK(a_data.size() == na * width);
  auto f = OMEGA_H_LAMBDA(LO a) {
    for (auto b = a2b[a]; b < a2b[a + 1]; ++b) {
      for (Int j = 0; j < width; ++j) {
        b_data[b * width + j] = a_data[a * width + j];
      }
    }
  };
  parallel_for(na, f, "expand_into");
}

template <typename T>
Read<T> expand(Read<T> a_data, LOs a2b, Int width) {
  auto nb = a2b.last();
  Write<T> b_data(nb * width);
  expand_into(a_data, a2b, b_data, width);
  return b_data;
}

template <typename T>
Read<T> permute(Read<T> a_data, LOs a2b, Int width) {
  auto nb = a2b.size();
  Write<T> b_data(nb * width);
  map_into(a_data, a2b, b_data, width);
  return b_data;
}

LOs multiply_fans(LOs a2b, LOs a2c) {
  auto b_degrees = get_degrees(a2b);
  auto c_degrees = get_degrees(a2c);
  LOs degrees = multiply_each(b_degrees, c_degrees);
  auto a2bc = offset_scan(degrees);
  return a2bc;
}

LOs compound_maps(LOs a2b, LOs b2c) {
  LO na = a2b.size();
  Write<LO> a2c(a2b.size());
  auto f = OMEGA_H_LAMBDA(LO a) {
    LO b = a2b[a];
    LO c = b2c[b];
    a2c[a] = c;
  };
  parallel_for(na, f, "compound_maps");
  return a2c;
}

LOs invert_permutation(LOs a2b) {
  Write<LO> b2a(a2b.size());
  auto f = OMEGA_H_LAMBDA(LO a) { b2a[a2b[a]] = a; };
  parallel_for(a2b.size(), f, "invert_permutation");
  return b2a;
}

Read<I8> invert_marks(Read<I8> marks) {
  Write<I8> out(marks.size());
  auto f = OMEGA_H_LAMBDA(LO i) { out[i] = !marks[i]; };
  parallel_for(out.size(), f, "invert_marks");
  return out;
}

LOs collect_marked(Read<I8> marks) {
  auto ntotal = marks.size();
  auto offsets = offset_scan(marks);
  auto nmarked = offsets.last();
  Write<LO> marked(nmarked);
  auto f = OMEGA_H_LAMBDA(LO i) {
    if (marks[i]) marked[offsets[i]] = i;
  };
  parallel_for(ntotal, f, "collect_marked");
  return marked;
}

Read<I8> mark_image(LOs a2b, LO nb) {
  auto na = a2b.size();
  Write<I8> out(nb, 0);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto b = a2b[a];
    out[b] = 1;
  };
  parallel_for(na, f, "mark_image");
  return out;
}

LOs invert_injective_map(LOs a2b, LO nb) {
  Write<LO> b2a(nb, -1);
  auto f = OMEGA_H_LAMBDA(LO a) { b2a[a2b[a]] = a; };
  parallel_for(a2b.size(), f, "invert_injective_map");
  return b2a;
}

LOs invert_funnel(LOs ab2a, LO na) {
  LO nab = ab2a.size();
  Write<LO> a2ab(na + 1, -1);
  a2ab.set(0, 0);
  auto f = OMEGA_H_LAMBDA(LO ab) {
    LO a_end = ab2a[ab];
    LO a_start = ab2a[ab + 1];
    if (a_end != a_start) {
      a2ab[a_end + 1] = ab + 1;
    }
  };
  parallel_for(nab - 1, f, "invert_funnel");
  if (nab) {
    LO a_end = ab2a.get(nab - 1);
    a2ab.set(a_end + 1, nab);
  }
  fill_right(a2ab);
  return a2ab;
}

Graph invert_map_by_sorting(LOs a2b, LO nb) {
  auto& ab2b = a2b;
  auto ba2ab = sort_by_keys(ab2b);
  auto ba2b = unmap(ba2ab, ab2b, 1);
  auto b2ba = invert_funnel(ba2b, nb);
  auto& ba2a = ba2ab;
  return Graph(b2ba, ba2a);
}

Graph invert_map_by_atomics(LOs a2b, LO nb, std::string const& b2ba_name,
    std::string const& ba2a_name) {
  auto na = a2b.size();
  Write<LO> degrees(nb, 0);
  auto count = OMEGA_H_LAMBDA(LO a) { atomic_increment(&degrees[a2b[a]]); };
  parallel_for(na, count, "invert_map_by_atomics(count)");
  auto b2ba = offset_scan(Read<LO>(degrees), b2ba_name);
  auto nba = b2ba.get(nb);
  Write<LO> write_ba2a(nba, ba2a_name);
  auto positions = Write<LO>(nb, 0);
  auto fill = OMEGA_H_LAMBDA(LO a) {
    auto b = a2b[a];
    auto first = b2ba[b];
    auto j = atomic_fetch_add<LO>(&positions[b], 1);
    write_ba2a[first + j] = a;
  };
  parallel_for(na, fill, "invert_map_by_atomics(fill");
  auto ba2a = LOs(write_ba2a);
  return Graph(b2ba, ba2a);
}

LOs get_degrees(LOs offsets, std::string const& name) {
  Write<LO> degrees(offsets.size() - 1, name);
  auto f = OMEGA_H_LAMBDA(LO i) { degrees[i] = offsets[i + 1] - offsets[i]; };
  parallel_for(degrees.size(), f, "get_degrees");
  return degrees;
}

LOs invert_fan(LOs a2b) {
  auto na = a2b.size() - 1;
  auto nb = a2b.last();
  Write<LO> b2a(nb, -1);
  auto f = OMEGA_H_LAMBDA(LO a) {
    if (a2b[a] != a2b[a + 1]) b2a[a2b[a]] = a;
  };
  parallel_for(na, f, "invert_fan");
  fill_right(b2a);
  return b2a;
}

Bytes mark_fan_preimage(LOs a2b) {
  OMEGA_H_CHECK(a2b.size() >= 1);
  auto out = Write<Byte>(a2b.size() - 1);
  auto f = OMEGA_H_LAMBDA(LO i) { out[i] = (a2b[i] != a2b[i + 1]); };
  parallel_for(out.size(), f);
  return out;
}

template <typename Functor>
static Read<typename Functor::input_type> fan_reduce_tmpl(
    LOs a2b, Read<typename Functor::input_type> b_data, Int width) {
  using T = typename Functor::input_type;
  using VT = typename Functor::value_type;
  OMEGA_H_CHECK(a2b.last() * width == b_data.size());
  auto na = a2b.size() - 1;
  Write<T> a_data(na * width);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto functor = Functor();
    for (Int j = 0; j < width; ++j) {
      VT res;
      functor.init(res);
      for (auto b = a2b[a]; b < a2b[a + 1]; ++b) {
        VT update = b_data[b * width + j];
        functor.join(res, update);
      }
      a_data[a * width + j] = static_cast<T>(res);
    }
  };
  parallel_for(na, f, "fan_reduce");
  return a_data;
}

template <typename T>
Read<T> fan_reduce(LOs a2b, Read<T> b_data, Int width, Omega_h_Op op) {
  switch (op) {
    case OMEGA_H_MIN:
      return fan_reduce_tmpl<MinFunctor<T>>(a2b, b_data, width);
    case OMEGA_H_MAX:
      return fan_reduce_tmpl<MaxFunctor<T>>(a2b, b_data, width);
    case OMEGA_H_SUM:
      return fan_reduce_tmpl<SumFunctor<T>>(a2b, b_data, width);
  }
  OMEGA_H_NORETURN(Read<T>());
}

#define INST_T(T)                                                              \
  template void add_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width); \
  template void map_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width); \
  template void map_into_range(                                                \
      Read<T> a_data, LO begin, LO end, Write<T> b_data, Int width);           \
  template Read<T> map_onto(Read<T> a_data, LOs a2b, LO nb, T, Int width);     \
  template Write<T> unmap(LOs a2b, Read<T> b_data, Int width);                 \
  template Read<T> unmap_range(LO begin, LO end, Read<T> b_data, Int width);   \
  template Read<T> expand(Read<T> a_data, LOs a2b, Int width);                 \
  template void expand_into(                                                   \
      Read<T> a_data, LOs a2b, Write<T> b_data, Int width);                    \
  template Read<T> permute(Read<T> a_data, LOs a2b, Int width);                \
  template Read<T> fan_reduce(                                                 \
      LOs a2b, Read<T> b_data, Int width, Omega_h_Op op);
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T

}  // end namespace Omega_h
