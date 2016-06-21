template <Int deg>
struct IsMatch;

template <>
struct IsMatch<2> {
  DEVICE static bool eval(
      LOs const& av2v,
      LO a_begin,
      LOs const& bv2v,
      LO b_begin,
      Int which_down,
      I8& match_code) {
    if (av2v[a_begin + 1] == bv2v[b_begin + (1 - which_down)]) {
      match_code = make_code(false, which_down, 0);
      return true;
    }
    return false;
  }
};

template <>
struct IsMatch<3> {
  DEVICE static bool eval(
      LOs const& av2v,
      LO a_begin,
      LOs const& bv2v,
      LO b_begin,
      Int which_down,
      I8& match_code) {
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 1) % 3)] &&
        av2v[a_begin + 2] == bv2v[b_begin + ((which_down + 2) % 3)]) {
      match_code = make_code(false, rotation_to_first<3>(which_down), 0);
      return true;
    }
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 2) % 3)] &&
        av2v[a_begin + 2] == bv2v[b_begin + ((which_down + 1) % 3)]) {
      match_code = make_code(true, rotation_to_first<3>(which_down), 0);
      return true;
    }
    return false;
  }
};

template <Int deg>
void find_matches(LOs av2v, LOs bv2v, Adj v2b,
    LOs* a2b_out, Read<I8>* codes_out) {
  LO na = av2v.size() / deg;
  LOs v2vb = v2b.a2ab;
  LOs vb2b = v2b.ab2b;
  Read<I8> vb_codes = v2b.codes;
  Write<LO> a2b(na);
  Write<I8> codes(na);
  auto f = LAMBDA(LO a) {
    LO a_begin = a * deg;
    LO v0 = av2v[a_begin];
    LO vb_begin = v2vb[v0];
    LO vb_end = v2vb[v0 + 1];
    for (LO vb = vb_begin; vb < vb_end; ++vb) {
      LO b = vb2b[vb];
      auto vb_code = vb_codes[vb];
      Int which_down = code_which_down(vb_code);
      LO b_begin = b * deg;
      I8 match_code;
      if (IsMatch<deg>::eval(av2v, a_begin, bv2v, b_begin,
            which_down, match_code)) {
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

void find_matches(Int dim, LOs av2v, LOs bv2v, Adj v2b,
    LOs* a2b_out, Read<I8>* codes_out) {
  if (dim == 1)
    find_matches<2>(av2v, bv2v, v2b, a2b_out, codes_out);
  if (dim == 2)
    find_matches<3>(av2v, bv2v, v2b, a2b_out, codes_out);
}

Adj reflect_down(LOs hv2v, LOs lv2v, Adj v2l,
    Int high_dim, Int low_dim) {
  LOs uv2v = form_uses(hv2v, high_dim, low_dim);
  LOs hl2l;
  Read<I8> codes;
  find_matches(low_dim, uv2v, lv2v, v2l, &hl2l, &codes);
  return Adj(hl2l, codes);
}

Adj reflect_down(LOs hv2v, LOs lv2v, LO nv,
    Int high_dim, Int low_dim) {
  Int nverts_per_low = simplex_degrees[low_dim][0];
  LO nl = lv2v.size() / nverts_per_low;
  auto l2v = Adj(lv2v);
  Adj v2l = invert(l2v, nverts_per_low, nv, Read<GO>(nl, 0, 1));
  return reflect_down(hv2v, lv2v, v2l, high_dim, low_dim);
}
