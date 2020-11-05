#include "Omega_h_refine_topology.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

void refine_edges_to_pairs(Mesh* mesh, LOs keys2edges, LOs keys2midverts,
    LOs old_verts2new_verts, LOs& keys2pairs, LOs& pair_verts2verts) {
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto nkeys = keys2edges.size();
  auto ndoms = nkeys;
  auto npairs = ndoms * 2;
  Write<LO> pair_verts2verts_w(npairs * 2);
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto edge = keys2edges[key];
    auto midvert = keys2midverts[key];
    pair_verts2verts_w[key * 4 + 0] =
        old_verts2new_verts[edge_verts2verts[edge * 2 + 0]];
    pair_verts2verts_w[key * 4 + 1] = midvert;
    pair_verts2verts_w[key * 4 + 2] = midvert;
    pair_verts2verts_w[key * 4 + 3] =
        old_verts2new_verts[edge_verts2verts[edge * 2 + 1]];
  };
  parallel_for(nkeys, f, "refine_edges_to_pairs");
  keys2pairs = LOs(nkeys + 1, 0, 2);
  pair_verts2verts = pair_verts2verts_w;
}

/* "domains" (doms for short) are the (dim)-dimensional
   entities adjacent to a key edge.
   each one of them is split into two "pair" entities
   and a "cut" entity running down the middle.
   when (dim == 1), the "cut" entities are vertices
   so some special handling is needed in that case
   (see refine_edge_interiors() above).

   both pairs and cuts are part of the general set of
   "product" entities (prods for short), i.e. new cavity
   interior entities created during mesh modification.
 */

void refine_domains_to_pairs(Mesh* mesh, Int dim, LOs keys2edges,
    LOs keys2midverts, LOs old_verts2new_verts, LOs& keys2pairs,
    LOs& pair_verts2verts) {
  OMEGA_H_CHECK(dim > VERT);
  if (dim == EDGE) {
    refine_edges_to_pairs(mesh, keys2edges, keys2midverts, old_verts2new_verts,
        keys2pairs, pair_verts2verts);
    return;
  }
  auto nkeys = keys2edges.size();
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto dom_verts2verts = mesh->ask_verts_of(dim);
  auto edges2doms = mesh->ask_up(EDGE, dim);
  auto edges2edge_doms = edges2doms.a2ab;
  auto edge_doms2doms = edges2doms.ab2b;
  auto edge_dom_codes = edges2doms.codes;
  auto edge_dom_degrees = get_degrees(edges2edge_doms);
  auto key_dom_degrees = read(unmap(keys2edges, edge_dom_degrees, 1));
  auto keys2key_doms = offset_scan(key_dom_degrees);
  auto ndoms = keys2key_doms.last();
  auto npairs = ndoms * 2;
  keys2pairs = multiply_each_by(keys2key_doms, 2);
  Write<LO> pair_verts2verts_w(npairs * (dim + 1));
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto edge = keys2edges[key];
    auto midvert = keys2midverts[key];
    auto pair = keys2pairs[key];
    for (auto edge_dom = edges2edge_doms[edge];
         edge_dom < edges2edge_doms[edge + 1]; ++edge_dom) {
      auto dom = edge_doms2doms[edge_dom];
      auto code = edge_dom_codes[edge_dom];
      auto dde = code_which_down(code);
      auto rot = code_rotation(code);
      for (Int eev = 0; eev < 2; ++eev) {
        /* a new cell is formed from an old cell by finding
           its side that is opposite to one of the edge endpoints
           and connecting it to the midpoint to form the new cell */
        auto dev = eev ^ rot;
        auto ddv = simplex_down_template(dim, EDGE, dde, dev);
        auto dds = simplex_opposite_template(dim, VERT, ddv);
        auto ppv2v = &pair_verts2verts_w[pair * (dim + 1)];
        for (Int dsv = 0; dsv < dim; ++dsv) {
          auto ddv2 = simplex_down_template(dim, dim - 1, dds, dsv);
          auto ov = dom_verts2verts[dom * (dim + 1) + ddv2];
          auto nv = old_verts2new_verts[ov];
          ppv2v[dsv] = nv;
        }
        ppv2v[dim] = midvert;
        flip_new_elem(dim, ppv2v);
        ++pair;
      }
    }
  };
  parallel_for(nkeys, f, "refine_domains_to_pairs");
  pair_verts2verts = pair_verts2verts_w;
}

void refine_domains_to_cuts(Mesh* mesh, Int dim, LOs keys2edges,
    LOs keys2midverts, LOs old_verts2new_verts, LOs& keys2cuts,
    LOs& cut_verts2verts) {
  OMEGA_H_CHECK(dim > EDGE);
  auto nkeys = keys2edges.size();
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto dom_verts2verts = mesh->ask_verts_of(dim);
  auto edges2doms = mesh->ask_up(EDGE, dim);
  auto edges2edge_doms = edges2doms.a2ab;
  auto edge_doms2doms = edges2doms.ab2b;
  auto edge_dom_codes = edges2doms.codes;
  auto edge_dom_degrees = get_degrees(edges2edge_doms);
  auto key_dom_degrees = read(unmap(keys2edges, edge_dom_degrees, 1));
  auto keys2key_doms = offset_scan(key_dom_degrees);
  auto ndoms = keys2key_doms.last();
  auto ncuts = ndoms;
  keys2cuts = keys2key_doms;
  Write<LO> cut_verts2verts_w(ncuts * (dim));
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto edge = keys2edges[key];
    auto midvert = keys2midverts[key];
    auto cut = keys2cuts[key];
    for (auto edge_dom = edges2edge_doms[edge];
         edge_dom < edges2edge_doms[edge + 1]; ++edge_dom) {
      auto dom = edge_doms2doms[edge_dom];
      auto code = edge_dom_codes[edge_dom];
      auto dde = code_which_down(code);
      /* a "cut" is formed by connecting its opposite
         "tip" (running out of words) to the new midpoint
         vertex. for triangle domains, the tip is the vertex
         not adjacent to the key edge. for tet domains, the tip
         is the edge not adjacent to the key edge. */
      auto ccv2v = &cut_verts2verts_w[cut * dim];
      auto ddt = simplex_opposite_template(dim, EDGE, dde);
      for (Int dtv = 0; dtv < dim - 1; ++dtv) {
        auto ddv2 = simplex_down_template(dim, dim - 2, ddt, dtv);
        auto ov = dom_verts2verts[dom * (dim + 1) + ddv2];
        auto nv = old_verts2new_verts[ov];
        ccv2v[dtv] = nv;
      }
      ccv2v[dim - 1] = midvert;
      ++cut;
    }
  };
  parallel_for(nkeys, f, "refine_domains_to_cuts");
  cut_verts2verts = cut_verts2verts_w;
}

void combine_pairs_and_cuts(Int ent_dim, LOs keys2cuts, LOs keys2pairs,
    LOs cut_verts2verts, LOs pair_verts2verts, LOs& keys2prods,
    LOs& prod_verts2verts) {
  auto nkeys = keys2cuts.size() - 1;
  OMEGA_H_CHECK(nkeys == keys2pairs.size() - 1);
  auto keys2ncuts = get_degrees(keys2cuts);
  auto keys2npairs = get_degrees(keys2pairs);
  auto keys2nprods = add_each(keys2ncuts, keys2npairs);
  keys2prods = offset_scan(keys2nprods);
  auto nprods = keys2prods.last();
  auto nppv = simplex_degree(ent_dim, VERT);
  auto prod_verts2verts_w = Write<LO>(nprods * nppv);
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto prod = keys2prods[key];
    for (auto pair = keys2pairs[key]; pair < keys2pairs[key + 1]; ++pair) {
      for (Int ppv = 0; ppv < nppv; ++ppv) {
        prod_verts2verts_w[prod * nppv + ppv] =
            pair_verts2verts[pair * nppv + ppv];
      }
      ++prod;
    }
    for (auto cut = keys2cuts[key]; cut < keys2cuts[key + 1]; ++cut) {
      for (Int ppv = 0; ppv < nppv; ++ppv) {
        prod_verts2verts_w[prod * nppv + ppv] =
            cut_verts2verts[cut * nppv + ppv];
      }
      ++prod;
    }
  };
  parallel_for(nkeys, f, "combine_pairs_and_cuts");
  prod_verts2verts = prod_verts2verts_w;
}

void refine_products(Mesh* mesh, Int ent_dim, LOs keys2edges, LOs keys2midverts,
    LOs old_verts2new_verts, LOs& keys2prods, LOs& prod_verts2verts) {
  auto keys2pairs = LOs();
  auto pair_verts2verts = LOs();
  refine_domains_to_pairs(mesh, ent_dim, keys2edges, keys2midverts,
      old_verts2new_verts, keys2pairs, pair_verts2verts);
  if (ent_dim == mesh->dim()) {
    keys2prods = keys2pairs;
    prod_verts2verts = pair_verts2verts;
  } else {
    auto keys2cuts = LOs();
    auto cut_verts2verts = LOs();
    refine_domains_to_cuts(mesh, ent_dim + 1, keys2edges, keys2midverts,
        old_verts2new_verts, keys2cuts, cut_verts2verts);
    combine_pairs_and_cuts(ent_dim, keys2cuts, keys2pairs, cut_verts2verts,
        pair_verts2verts, keys2prods, prod_verts2verts);
  }
}

}  // end namespace Omega_h
