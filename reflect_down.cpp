#include "adjacency.hpp"

#include "align.hpp"
#include "array.hpp"
#include "loop.hpp"
#include "simplices.hpp"

namespace osh {

template <Int deg>
struct IsMatch;

template <>
struct IsMatch<2> {
  template <typename T>
  DEVICE static bool eval(Read<T> const& av2v, LO a_begin, Read<T> const& bv2v,
                          LO b_begin, Int which_down, I8* match_code) {
    if (av2v[a_begin + 1] == bv2v[b_begin + (1 - which_down)]) {
      *match_code = make_code(false, which_down, 0);
      return true;
    }
    return false;
  }
};

template <>
struct IsMatch<3> {
  template <typename T>
  DEVICE static bool eval(Read<T> const& av2v, LO a_begin, Read<T> const& bv2v,
                          LO b_begin, Int which_down, I8* match_code) {
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 1) % 3)] &&
        av2v[a_begin + 2] == bv2v[b_begin + ((which_down + 2) % 3)]) {
      *match_code = make_code(false, rotation_to_first<3>(which_down), 0);
      return true;
    }
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 2) % 3)] &&
        av2v[a_begin + 2] == bv2v[b_begin + ((which_down + 1) % 3)]) {
      *match_code = make_code(true, rotation_to_first<3>(which_down), 0);
      return true;
    }
    return false;
  }
};

template <Int deg, typename T>
void find_matches_deg(LOs a2fv, Read<T> av2v, Read<T> bv2v, Adj v2b,
                      LOs* a2b_out, Read<I8>* codes_out) {
  LO na = a2fv.size();
  CHECK(na * deg == av2v.size());
  LOs v2vb = v2b.a2ab;
  LOs vb2b = v2b.ab2b;
  Read<I8> vb_codes = v2b.codes;
  Write<LO> a2b(na);
  Write<I8> codes(na);
  auto f = LAMBDA(LO a) {
    auto fv = a2fv[a];
    auto a_begin = a * deg;
    auto vb_begin = v2vb[fv];
    auto vb_end = v2vb[fv + 1];
    for (LO vb = vb_begin; vb < vb_end; ++vb) {
      auto b = vb2b[vb];
      auto vb_code = vb_codes[vb];
      auto which_down = code_which_down(vb_code);
      auto b_begin = b * deg;
      I8 match_code;
      if (IsMatch<deg>::eval(av2v, a_begin, bv2v, b_begin, which_down,
                             &match_code)) {
        a2b[a] = b;
        codes[a] = match_code;
        return;
      }
    }
    NORETURN();
  };
  parallel_for(na, f);
  *a2b_out = a2b;
  *codes_out = codes;
}

template <typename T>
void find_matches_ex(Int deg, LOs a2fv, Read<T> av2v, Read<T> bv2v, Adj v2b,
                     LOs* a2b_out, Read<I8>* codes_out) {
  if (deg == 2) {
    find_matches_deg<2>(a2fv, av2v, bv2v, v2b, a2b_out, codes_out);
  } else if (deg == 3) {
    find_matches_deg<3>(a2fv, av2v, bv2v, v2b, a2b_out, codes_out);
  }
}

void find_matches(Int dim, LOs av2v, LOs bv2v, Adj v2b, LOs* a2b_out,
                  Read<I8>* codes_out) {
  auto deg = dim + 1;
  auto a2fv = get_component(av2v, deg, 0);
  find_matches_ex(deg, a2fv, av2v, bv2v, v2b, a2b_out, codes_out);
}

Adj reflect_down(LOs hv2v, LOs lv2v, Adj v2l, Int high_dim, Int low_dim) {
  LOs uv2v = form_uses(hv2v, high_dim, low_dim);
  LOs hl2l;
  Read<I8> codes;
  find_matches(low_dim, uv2v, lv2v, v2l, &hl2l, &codes);
  return Adj(hl2l, codes);
}

Adj reflect_down(LOs hv2v, LOs lv2v, LO nv, Int high_dim, Int low_dim) {
  Int nverts_per_low = simplex_degrees[low_dim][0];
  LO nl = lv2v.size() / nverts_per_low;
  auto l2v = Adj(lv2v);
  Adj v2l = invert(l2v, nverts_per_low, nv, Read<GO>(nl, 0, 1));
  return reflect_down(hv2v, lv2v, v2l, high_dim, low_dim);
}

#define INST(T) \
template \
void find_matches_ex(Int deg, LOs a2fv, Read<T> av2v, Read<T> bv2v, Adj v2b, \
                     LOs* a2b_out, Read<I8>* codes_out);
INST(LO)
INST(GO)
#undef INST

} //end namespace osh
