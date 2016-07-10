#include "internal.hpp"

namespace osh {

template <Int deg>
struct IsFlipped;
template <>
struct IsFlipped<3> {
  template <typename T>
  INLINE static bool is(T adj[]) {
    return adj[2] < adj[1];
  }
};
template <>
struct IsFlipped<2> {
  template <typename T>
  INLINE static bool is(T adj[]) {
    (void)adj;
    return false;
  }
};

template <Int deg, typename T>
static Read<I8> get_codes_to_canonical_deg(Read<T> ev2v) {
  LO nev = ev2v.size();
  LO ne = nev / deg;
  Write<I8> codes(ne);
  auto f = LAMBDA(LO e) {
    LO begin = e * deg;
    /* find the smallest vertex */
    Int min_j = 0;
    auto min_v = ev2v[begin];
    for (Int j = 1; j < deg; ++j) {
      LO ev = j + begin;
      auto v = ev2v[ev];
      if (v < min_v) {
        min_j = j;
        min_v = v;
      }
    }
    /* rotate to make it first */
    auto rotation = rotation_to_first<deg>(min_j);
    T tmp[deg];
    rotate_adj<deg>(rotation, &ev2v[begin], tmp);
    auto is_flipped = IsFlipped<deg>::is(tmp);
    codes[e] = make_code(is_flipped, rotation, 0);
  };
  parallel_for(ne, f);
  return codes;
}

template <typename T>
Read<I8> get_codes_to_canonical(Int deg, Read<T> ev2v) {
  if (deg == 3) {
    return get_codes_to_canonical_deg<3>(ev2v);
  } else {
    CHECK(deg == 2);
    return get_codes_to_canonical_deg<2>(ev2v);
  }
}

#define INST(T) template Read<I8> get_codes_to_canonical(Int deg, Read<T> ev2v);
INST(LO)
INST(GO)
#undef INST

/* check whether adjacent lists of (deg) vertices
   are the same */
DEVICE static bool are_equal(Int deg, LOs const& canon, LO e0, LO e1) {
  LO a = e0 * deg;
  LO b = e1 * deg;
  for (LO j = 0; j < deg; ++j)
    if (canon[a + j] != canon[b + j]) return false;
  return true;
}

Read<I8> find_canonical_jumps(Int deg, LOs canon, LOs e_sorted2e) {
  auto ne = e_sorted2e.size();
  Write<I8> jumps(ne, 0);
  auto f = LAMBDA(LO e_sorted) {
    auto e0 = e_sorted2e[e_sorted];
    auto e1 = e_sorted2e[e_sorted + 1];
    if (!are_equal(deg, canon, e0, e1)) jumps[e_sorted] = 1;
  };
  parallel_for(ne - 1, f);
  if (jumps.size()) jumps.set(jumps.size() - 1, 1);
  return jumps;
}

static LOs find_unique_deg(Int deg, LOs uv2v) {
  auto codes = get_codes_to_canonical(deg, uv2v);
  auto uv2v_canon = align_ev2v(deg, uv2v, codes);
  auto sorted2u = sort_by_keys(uv2v_canon, deg);
  auto jumps = find_canonical_jumps(deg, uv2v_canon, sorted2u);
  auto e2sorted = collect_marked(jumps);
  auto e2u = compound_maps(e2sorted, sorted2u);
  return unmap<LO>(e2u, uv2v, deg);
}

LOs find_unique(LOs hv2v, Int high_dim, Int low_dim) {
  auto uv2v = form_uses(hv2v, high_dim, low_dim);
  return find_unique_deg(low_dim + 1, uv2v);
}

} //end namespace osh
