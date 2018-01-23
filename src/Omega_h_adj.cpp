#include "Omega_h_adj.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_control.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_scan.hpp"
#include "Omega_h_simplex.hpp"
#include "Omega_h_sort.hpp"
#include "Omega_h_timer.hpp"

#include <iostream>

namespace Omega_h {

Adj unmap_adjacency(LOs a2b, Adj b2c) {
  auto b2bc = b2c.a2ab;
  auto bc2c = b2c.ab2b;
  auto bc_codes = b2c.codes;
  auto b_degrees = get_degrees(b2bc);
  auto a_degrees = unmap(a2b, b_degrees, 1);
  auto a2ac = offset_scan(a_degrees);
  auto na = a2b.size();
  auto nac = a2ac.last();
  Write<LO> ac2c(nac);
  auto ac_codes = Write<I8>(nac);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto b = a2b[a];
    auto bc = b2bc[b];
    for (auto ac = a2ac[a]; ac < a2ac[a + 1]; ++ac) {
      ac2c[ac] = bc2c[bc];
      ac_codes[ac] = bc_codes[bc];
      ++bc;
    }
  };
  parallel_for(na, f, "unmap_adjacency");
  return Adj(a2ac, ac2c, ac_codes);
}

template <Int deg>
struct IsFlipped;
template <>
struct IsFlipped<3> {
  template <typename T>
  OMEGA_H_INLINE static bool is(T adj[]) {
    return adj[2] < adj[1];
  }
};

template <>
struct IsFlipped<2> {
  template <typename T>
  OMEGA_H_INLINE static bool is(T adj[]) {
    (void)adj;
    return false;
  }
};

template <Int deg, typename T>
static Read<I8> get_codes_to_canonical_deg(Read<T> ev2v) {
  auto nev = ev2v.size();
  auto ne = nev / deg;
  Write<I8> codes(ne);
  auto f = OMEGA_H_LAMBDA(LO e) {
    auto begin = e * deg;
    /* find the smallest vertex */
    Int min_j = 0;
    auto min_v = ev2v[begin];
    for (Int j = 1; j < deg; ++j) {
      auto ev = j + begin;
      auto v = ev2v[ev];
      if (v < min_v) {
        min_j = j;
        min_v = v;
      }
    }
    /* rotate to make it first */
    auto rotation = rotation_to_first<deg>(min_j);
    T tmp[deg];
    rotate_adj<deg>(rotation, ev2v, begin, tmp, 0);
    auto is_flipped = IsFlipped<deg>::is(tmp);
    codes[e] = make_code(is_flipped, rotation, 0);
  };
  parallel_for(ne, f, "get_codes_to_canonical");
  return codes;
}

template <typename T>
Read<I8> get_codes_to_canonical(Int deg, Read<T> ev2v) {
  if (deg == 3) return get_codes_to_canonical_deg<3>(ev2v);
  if (deg == 2) return get_codes_to_canonical_deg<2>(ev2v);
  OMEGA_H_NORETURN(Read<I8>());
}

/* check whether adjacent lists of (deg) vertices
   are the same */
OMEGA_H_DEVICE static bool are_equal(Int deg, LOs const& canon, LO e0, LO e1) {
  auto a = e0 * deg;
  auto b = e1 * deg;
  for (LO j = 0; j < deg; ++j) {
    if (canon[a + j] != canon[b + j]) return false;
  }
  return true;
}

Read<I8> find_canonical_jumps(Int deg, LOs canon, LOs e_sorted2e) {
  auto ne = e_sorted2e.size();
  Write<I8> jumps(ne, 0);
  auto f = OMEGA_H_LAMBDA(LO e_sorted) {
    auto e0 = e_sorted2e[e_sorted];
    auto e1 = e_sorted2e[e_sorted + 1];
    if (!are_equal(deg, canon, e0, e1)) jumps[e_sorted] = 1;
  };
  parallel_for(ne - 1, f, "find_canonical_jumps");
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

LOs form_uses(LOs hv2v, Int high_dim, Int low_dim) {
  Int nverts_per_high = simplex_degree(high_dim, 0);
  Int nverts_per_low = simplex_degree(low_dim, 0);
  Int nlows_per_high = simplex_degree(high_dim, low_dim);
  LO nhigh = hv2v.size() / nverts_per_high;
  LO nuses = nhigh * nlows_per_high;
  Write<LO> uv2v(nuses * nverts_per_low);
  auto f = OMEGA_H_LAMBDA(LO h) {
    LO h_begin = h * nverts_per_high;
    for (Int u = 0; u < nlows_per_high; ++u) {
      LO u_begin = (h * nlows_per_high + u) * nverts_per_low;
      for (Int uv = 0; uv < nverts_per_low; ++uv) {
        uv2v[u_begin + uv] =
            hv2v[h_begin + simplex_down_template(high_dim, low_dim, u, uv)];
      }
    }
  };
  parallel_for(nhigh, f, "form_uses");
  return uv2v;
}

static void sort_by_high_index(LOs l2lh, Write<LO> lh2h, Write<I8> codes) {
  LO nl = l2lh.size() - 1;
  auto f = OMEGA_H_LAMBDA(LO l) {
    LO begin = l2lh[l];
    LO end = l2lh[l + 1];
    for (LO j = begin; j < end; ++j) {
      LO k_min = j;
      GO min_h = lh2h[j];
      for (LO k = j + 1; k < end; ++k) {
        GO h = lh2h[k];
        if (h < min_h) {
          k_min = k;
          min_h = h;
        }
      }
      swap2(lh2h[j], lh2h[k_min]);
      swap2(codes[j], codes[k_min]);
    }
  };
  parallel_for(nl, f, "sort_by_high_index");
}

Adj invert_adj(Adj down, Int nlows_per_high, LO nlows) {
  begin_code("invert_adj");
  auto l2hl = invert_map_by_atomics(down.ab2b, nlows);
  auto l2lh = l2hl.a2ab;
  auto lh2hl = l2hl.ab2b;
  LO nlh = lh2hl.size();
  Read<I8> down_codes(down.codes);
  Write<LO> lh2h(nlh);
  Write<I8> codes(nlh);
  if (down_codes.exists()) {
    auto f = OMEGA_H_LAMBDA(LO lh) {
      LO hl = lh2hl[lh];
      LO h = hl / nlows_per_high;
      lh2h[lh] = h;
      Int which_down = hl % nlows_per_high;
      auto down_code = down_codes[hl];
      bool is_flipped = code_is_flipped(down_code);
      Int rotation = code_rotation(down_code);
      codes[lh] = make_code(is_flipped, rotation, which_down);
    };
    parallel_for(nlh, f, "full_codes");
  } else {
    auto f = OMEGA_H_LAMBDA(LO lh) {
      LO hl = lh2hl[lh];
      LO h = hl / nlows_per_high;
      lh2h[lh] = h;
      Int which_down = hl % nlows_per_high;
      codes[lh] = make_code(false, 0, which_down);
    };
    parallel_for(nlh, f, "easy_codes");
  }
  sort_by_high_index(l2lh, lh2h, codes);
  end_code();
  return Adj(l2lh, lh2h, codes);
}

template <Int deg>
struct IsMatch;

template <>
struct IsMatch<2> {
  template <typename T>
  OMEGA_H_DEVICE static bool eval(Read<T> const& av2v, LO a_begin,
      Read<T> const& bv2v, LO b_begin, Int which_down, I8* match_code) {
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
  OMEGA_H_DEVICE static bool eval(Read<T> const& av2v, LO a_begin,
      Read<T> const& bv2v, LO b_begin, Int which_down, I8* match_code) {
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

template <>
struct IsMatch<4> {
  template <typename T>
  OMEGA_H_DEVICE static bool eval(Read<T> const& av2v, LO a_begin,
      Read<T> const& bv2v, LO b_begin, Int which_down, I8* match_code) {
    if (av2v[a_begin + 2] != bv2v[b_begin + 2]) return false;
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 1) % 4)] &&
        av2v[a_begin + 3] == bv2v[b_begin + ((which_down + 3) % 4)]) {
      *match_code = 0;
      return true;
    }
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 3) % 4)] &&
        av2v[a_begin + 3] == bv2v[b_begin + ((which_down + 1) % 4)]) {
      *match_code = 0;
      return true;
    }
    return false;
  }
};

template <Int deg, typename T>
static void find_matches_deg(LOs a2fv, Read<T> av2v, Read<T> bv2v, Adj v2b,
    LOs* a2b_out, Read<I8>* codes_out, bool allow_duplicates) {
  LO na = a2fv.size();
  OMEGA_H_CHECK(na * deg == av2v.size());
  LOs v2vb = v2b.a2ab;
  LOs vb2b = v2b.ab2b;
  Read<I8> vb_codes = v2b.codes;
  Write<LO> a2b(na);
  Write<I8> codes(na);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto fv = a2fv[a];
    auto a_begin = a * deg;
    auto vb_begin = v2vb[fv];
    auto vb_end = v2vb[fv + 1];
    bool found = false;
    for (LO vb = vb_begin; vb < vb_end; ++vb) {
      auto b = vb2b[vb];
      auto vb_code = vb_codes[vb];
      auto which_down = code_which_down(vb_code);
      auto b_begin = b * deg;
      I8 match_code;
      if (IsMatch<deg>::eval(
              av2v, a_begin, bv2v, b_begin, which_down, &match_code)) {
        OMEGA_H_CHECK(!found);  // there can't be more than one!
        a2b[a] = b;
        codes[a] = match_code;
        found = true;
        if (allow_duplicates) break;
      }
    }
    OMEGA_H_CHECK(found);  // there can't be less than one!
  };
  parallel_for(na, f, "find_matches");
  *a2b_out = a2b;
  *codes_out = codes;
}

template <typename T>
void find_matches_ex(Int deg, LOs a2fv, Read<T> av2v, Read<T> bv2v, Adj v2b,
    LOs* a2b_out, Read<I8>* codes_out, bool allow_duplicates) {
  if (deg == 2) {
    find_matches_deg<2>(
        a2fv, av2v, bv2v, v2b, a2b_out, codes_out, allow_duplicates);
  } else if (deg == 3) {
    find_matches_deg<3>(
        a2fv, av2v, bv2v, v2b, a2b_out, codes_out, allow_duplicates);
  } else if (deg == 4) {
    find_matches_deg<4>(
        a2fv, av2v, bv2v, v2b, a2b_out, codes_out, allow_duplicates);
  }
}

void find_matches(
    Int dim, LOs av2v, LOs bv2v, Adj v2b, LOs* a2b_out, Read<I8>* codes_out) {
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
  Int nverts_per_low = simplex_degree(low_dim, 0);
  auto l2v = Adj(lv2v);
  Adj v2l = invert_adj(l2v, nverts_per_low, nv);
  return reflect_down(hv2v, lv2v, v2l, high_dim, low_dim);
}

Adj transit(Adj h2m, Adj m2l, Int high_dim, Int low_dim) {
  OMEGA_H_CHECK(3 >= high_dim);
  auto mid_dim = low_dim + 1;
  OMEGA_H_CHECK(high_dim > mid_dim);
  OMEGA_H_CHECK(low_dim == 1 || low_dim == 0);
  auto hm2m = h2m.ab2b;
  auto m2hm_codes = h2m.codes;
  auto ml2l = m2l.ab2b;
  auto ml_codes = m2l.codes;
  auto nmids_per_high = simplex_degree(high_dim, mid_dim);
  auto nlows_per_mid = simplex_degree(mid_dim, low_dim);
  auto nlows_per_high = simplex_degree(high_dim, low_dim);
  auto nhighs = hm2m.size() / nmids_per_high;
  Write<LO> hl2l(nhighs * nlows_per_high);
  Write<I8> codes;
  if (low_dim == 1) codes = Write<I8>(hl2l.size());
  auto f = OMEGA_H_LAMBDA(LO h) {
    auto hl_begin = h * nlows_per_high;
    auto hm_begin = h * nmids_per_high;
    for (Int hl = 0; hl < nlows_per_high; ++hl) {
      auto ut = simplex_up_template(high_dim, low_dim, hl, 0);
      auto hm = ut.up;
      auto hml = ut.which_down;
      auto m = hm2m[hm_begin + hm];
      auto m2hm_code = m2hm_codes[hm_begin + hm];
      auto hm2m_code = invert_alignment(nlows_per_mid, m2hm_code);
      auto ml = align_index(nlows_per_mid, low_dim, hml, hm2m_code);
      auto ml_begin = m * nlows_per_mid;
      auto l = ml2l[ml_begin + ml];
      // safety check for duplicates.
      // remove after this code is heavily exercised (or don't)
      for (Int hhl2 = 0; hhl2 < hl; ++hhl2) {
        OMEGA_H_CHECK(l != hl2l[hl_begin + hhl2]);
      }
      hl2l[hl_begin + hl] = l;
      if (low_dim == 1) {
        auto tet_tri_code = hm2m_code;
        auto tri_edge_code = ml_codes[ml_begin + ml];
        auto tet_tri_flipped = code_is_flipped(tet_tri_code);
        auto tri_edge_flipped = bool(code_rotation(tri_edge_code) == 1);
        auto canon_flipped = ut.is_flipped;
        bool tet_edge_flipped =
            tet_tri_flipped ^ tri_edge_flipped ^ canon_flipped;
        codes[hl_begin + hl] = make_code(false, tet_edge_flipped, 0);
      }
    }
  };
  parallel_for(nhighs, f, "transit");
  if (low_dim == 1) return Adj(hl2l, codes);
  return Adj(hl2l);
}

Graph verts_across_edges(Adj e2v, Adj v2e) {
  auto ev2v = e2v.ab2b;
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto v2ve_codes = v2e.codes;
  auto& v2vv = v2ve;
  Write<LO> vv2v(ve2e.size());
  auto f = OMEGA_H_LAMBDA(LO ve) {
    auto vv = ve;
    auto e = ve2e[ve];
    auto v2ve_code = v2ve_codes[ve];
    auto eev = code_which_down(v2ve_code);
    auto v = ev2v[e * 2 + (1 - eev)];
    vv2v[vv] = v;
  };
  parallel_for(vv2v.size(), f, "verts_across_edges");
  return Adj(v2vv, vv2v);
}

Graph edges_across_tris(Adj f2e, Adj e2f) {
  auto fe2e = f2e.ab2b;
  auto e2ef = e2f.a2ab;
  auto ef2f = e2f.ab2b;
  auto e2ef_codes = e2f.codes;
  auto ne = e2ef.size() - 1;
  auto e2ef_degrees = get_degrees(e2ef);
  auto e2ee_degrees = multiply_each_by(e2ef_degrees, 2);
  auto e2ee = offset_scan(e2ee_degrees);
  auto nee = e2ee.last();
  Write<LO> ee2e(nee);
  auto lambda = OMEGA_H_LAMBDA(LO e) {
    auto ef_begin = e2ef[e];
    auto ef_end = e2ef[e + 1];
    auto neef = ef_end - ef_begin;
    auto ee_begin = e2ee[e];
    for (Int eef = 0; eef < neef; ++eef) {
      auto ef = ef_begin + eef;
      auto f = ef2f[ef];
      auto e2ef_code = e2ef_codes[ef];
      auto ffe = code_which_down(e2ef_code);
      auto e1 = fe2e[f * 3 + ((ffe + 1) % 3)];
      auto e2 = fe2e[f * 3 + ((ffe + 2) % 3)];
      ee2e[ee_begin + eef * 2 + 0] = e1;
      ee2e[ee_begin + eef * 2 + 1] = e2;
    }
  };
  parallel_for(ne, lambda, "edges_across_tris");
  return Adj(e2ee, ee2e);
}

Graph edges_across_tets(Adj r2e, Adj e2r) {
  auto re2e = r2e.ab2b;
  auto e2er = e2r.a2ab;
  auto er2r = e2r.ab2b;
  auto e2er_codes = e2r.codes;
  auto ne = e2er.size() - 1;
  auto& e2ee = e2er;
  Write<LO> ee2e(er2r.size());
  auto f = OMEGA_H_LAMBDA(LO e) {
    auto er_begin = e2er[e];
    auto er_end = e2er[e + 1];
    for (auto er = er_begin; er < er_end; ++er) {
      auto r = er2r[er];
      auto e2er_code = e2er_codes[er];
      auto rre = code_which_down(e2er_code);
      auto rre_opp = simplex_opposite_template(REGION, EDGE, rre);
      auto re_begin = r * 6;
      auto e_opp = re2e[re_begin + rre_opp];
      auto ee = er;
      ee2e[ee] = e_opp;
    }
  };
  parallel_for(ne, f, "edges_across_tets");
  return Adj(e2ee, ee2e);
}

Graph elements_across_sides(
    Int dim, Adj elems2sides, Adj sides2elems, Read<I8> side_is_exposed) {
  auto elem_side2side = elems2sides.ab2b;
  auto side2side_elems = sides2elems.a2ab;
  auto side_elem2elem = sides2elems.ab2b;
  Int nsides_per_elem = dim + 1;
  auto nelems = elem_side2side.size() / nsides_per_elem;
  Write<LO> degrees(nelems);
  auto count = OMEGA_H_LAMBDA(LO elem) {
    auto begin = elem * nsides_per_elem;
    auto end = begin + nsides_per_elem;
    Int n = 0;
    for (auto elem_side = begin; elem_side < end; ++elem_side) {
      auto side = elem_side2side[elem_side];
      if (!side_is_exposed[side]) {
        OMEGA_H_CHECK(side2side_elems[side + 1] - side2side_elems[side] == 2);
        ++n;
      }
    }
    degrees[elem] = n;
  };
  parallel_for(nelems, count, "elements_across_sides(count)");
  auto elem2elem_elems = offset_scan(LOs(degrees));
  auto nelem_elems = elem2elem_elems.last();
  Write<LO> elem_elem2elem(nelem_elems);
  auto fill = OMEGA_H_LAMBDA(LO elem) {
    auto begin = elem * nsides_per_elem;
    auto end = begin + nsides_per_elem;
    LO elem_elem = elem2elem_elems[elem];
    for (auto elem_side = begin; elem_side < end; ++elem_side) {
      auto side = elem_side2side[elem_side];
      if (!side_is_exposed[side]) {
        auto side_elem = side2side_elems[side];
        if (side_elem2elem[side_elem] == elem) ++side_elem;
        elem_elem2elem[elem_elem] = side_elem2elem[side_elem];
        ++elem_elem;
      }
    }
  };
  parallel_for(nelems, fill, "elements_across_sides(fill)");
  return Graph(elem2elem_elems, elem_elem2elem);
}

#define INST(T)                                                                \
  template Read<I8> get_codes_to_canonical(Int deg, Read<T> ev2v);             \
  template void find_matches_ex(Int deg, LOs a2fv, Read<T> av2v, Read<T> bv2v, \
      Adj v2b, LOs* a2b_out, Read<I8>* codes_out, bool);
INST(LO)
INST(GO)
#undef INST

}  // namespace Omega_h
