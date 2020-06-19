#include "Omega_h_coarsen.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_collapse.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_simplex.hpp"

#include <iostream>

namespace Omega_h {

LOs get_verts_onto(Mesh* mesh, LOs rails2edges, Read<I8> rail_col_dirs) {
  auto nkeys = rails2edges.size();
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto keys2verts_onto_w = Write<LO>(nkeys, -1);
  auto set_key_onto = OMEGA_H_LAMBDA(LO key) {
    auto e = rails2edges[key];
    auto eev = rail_col_dirs[key];
    keys2verts_onto_w[key] = ev2v[e * 2 + (1 - eev)];
  };
  parallel_for(nkeys, set_key_onto, "get_verts_onto");
  return keys2verts_onto_w;
}

void mark_dead_ents(Mesh* mesh, LOs rails2edges, Read<I8> rail_col_dirs,
    Int cell_dim, Write<I8>& dead_cells, Write<I8>& dead_sides) {
  auto e2c = mesh->ask_up(EDGE, cell_dim);
  auto e2ec = e2c.a2ab;
  auto ec2c = e2c.ab2b;
  auto ec_codes = e2c.codes;
  auto cs2s = mesh->ask_down(cell_dim, cell_dim - 1).ab2b;
  auto nccs = simplex_degree(cell_dim, cell_dim - 1);
  auto nrails = rails2edges.size();
  auto f = OMEGA_H_LAMBDA(LO rail) {
    auto e = rails2edges[rail];
    auto eev_col = rail_col_dirs[rail];
    auto eev_onto = 1 - eev_col;
    for (auto ec = e2ec[e]; ec < e2ec[e + 1]; ++ec) {
      auto c = ec2c[ec];
      auto ec_code = ec_codes[ec];
      auto cce = code_which_down(ec_code);
      auto rot = code_rotation(ec_code);
      auto cev_onto = rot ^ eev_onto;
      auto ccv_onto = simplex_down_template(cell_dim, EDGE, cce, cev_onto);
      auto ccs_opp = simplex_opposite_template(cell_dim, VERT, ccv_onto);
      auto s_opp = cs2s[c * nccs + ccs_opp];
      dead_cells[c] = 1;
      dead_sides[s_opp] = 1;
    }
  };
  parallel_for(nrails, f, "mark_dead_ents");
}

HostFew<Read<I8>, 4> mark_dead_ents(
    Mesh* mesh, LOs rails2edges, Read<I8> rail_col_dirs) {
  HostFew<Write<I8>, 4> writes;
  writes[EDGE] = deep_copy(mark_image(rails2edges, mesh->nedges()));
  for (Int dim = EDGE + 1; dim <= mesh->dim(); ++dim)
    writes[dim] = Write<I8>(mesh->nents(dim), 0);
  for (Int dim = mesh->dim(); dim > EDGE; --dim)
    mark_dead_ents(
        mesh, rails2edges, rail_col_dirs, dim, writes[dim], writes[dim - 1]);
  HostFew<Read<I8>, 4> reads;
  for (Int dim = 0; dim < 4; ++dim) reads[dim] = writes[dim];
  return reads;
}

Adj find_coarsen_domains(
    Mesh* mesh, LOs keys2verts, Int ent_dim, Read<I8> ents_are_dead) {
  auto nkeys = keys2verts.size();
  auto v2e = mesh->ask_up(VERT, ent_dim);
  auto k2e = unmap_adjacency(keys2verts, v2e);
  auto k2ke = k2e.a2ab;
  auto ke2e = k2e.ab2b;
  auto ke_codes = k2e.codes;
  auto ke2k = invert_fan(k2ke);
  auto ents_are_live = invert_marks(ents_are_dead);
  auto kes_are_live = unmap(ke2e, ents_are_live, 1);
  auto lke2ke = collect_marked(kes_are_live);
  auto lke2k = unmap(lke2ke, ke2k, 1);
  auto lke_codes = unmap(lke2ke, ke_codes, 1);
  auto lke2e = unmap(lke2ke, ke2e, 1);
  auto k2lke = invert_funnel(lke2k, nkeys);
  return Adj(k2lke, lke2e, lke_codes);
}

LOs coarsen_topology(Mesh* mesh, LOs keys2verts_onto, Int dom_dim,
    Adj keys2doms, LOs old_verts2new_verts) {
  auto nccv = simplex_degree(dom_dim, VERT);
  auto cv2v = mesh->ask_verts_of(dom_dim);
  auto k2kc = keys2doms.a2ab;
  auto kc2c = keys2doms.ab2b;
  auto kc_codes = keys2doms.codes;
  auto ndoms = k2kc.last();
  auto nprods = ndoms;
  auto prod_verts2verts = Write<LO>(nprods * nccv);
  auto nkeys = keys2verts_onto.size();
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto v_onto = keys2verts_onto[key];
    for (auto kc = k2kc[key]; kc < k2kc[key + 1]; ++kc) {
      auto prod = kc;
      auto c = kc2c[kc];
      auto kc_code = kc_codes[kc];
      auto ccv_col = code_which_down(kc_code);
      auto ppv2v = &prod_verts2verts[prod * nccv];
      for (Int ccv = 0; ccv < nccv; ++ccv) {
        LO old_v;
        if (ccv == ccv_col) {
          old_v = v_onto;
        } else {
          old_v = cv2v[c * nccv + ccv];
        }
        ppv2v[ccv] = old_verts2new_verts[old_v];
      }
    }
  };
  parallel_for(nkeys, f, "coarsen_topology");
  return prod_verts2verts;
}

}  // end namespace Omega_h
