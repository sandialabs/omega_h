#include "Omega_h_collapse.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

Read<I8> check_collapse_class(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  auto ncands = cands2edges.size();
  auto verts2class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2class_dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto cand_codes_w = Write<I8>(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto code = cand_codes[cand];
    auto edge = cands2edges[cand];
    LO eev2v[2];
    for (Int eev = 0; eev < 2; ++eev)
      eev2v[eev] = edge_verts2verts[edge * 2 + eev];
    Int eev_cd[2];
    for (Int eev = 0; eev < 2; ++eev) eev_cd[eev] = verts2class_dim[eev2v[eev]];
    Int e_cd = edges2class_dim[edge];
    /* if both vertices are classified onto the same dimension,
       the edge must also be classified onto that dimension */
    if (eev_cd[0] == eev_cd[1]) {
      if (eev_cd[0] != e_cd) code = DONT_COLLAPSE;
    } else {
      for (Int col_v = 0; col_v < 2; ++col_v) {
        if (collapses(cand_codes[cand], col_v)) {
          /* otherwise, the vertex to collapse and the edge must
             be classified onto the same dimension */
          if (eev_cd[col_v] != e_cd) {
            code = dont_collapse(code, col_v);
          }
        }
      }
    }
    cand_codes_w[cand] = code;
  };
  parallel_for(ncands, f, "check_collapse_class");
  return cand_codes_w;
}

/* we also have to check that for every entity being
   collapsed, it is not exposing (re-classifying)
   any boundary entities.
   this is part of the overall axiom that cavity
   boundaries are preserved.

If we are
collapsing
this vertex ~~~~> *
                  | \
                  |   \~~~~~~~~~ This edge ("side" opposite the "onto" vertex)
                | |     \        must have the same classification as...
     Along this | |       \
     direction  | |         *
                V |       /
                  |  X~~~~~~~~~~ ...this triangle ("cell").
                  |   /
Onto this         | /
vertex ~~~~~~~~>  *

*/

Read<I8> check_collapse_exposure(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes, Int cell_dim) {
  auto e2c = mesh->ask_up(EDGE, cell_dim);
  auto e2ec = e2c.a2ab;
  auto ec2c = e2c.ab2b;
  auto ec_codes = e2c.codes;
  auto cs2s = mesh->ask_down(cell_dim, cell_dim - 1).ab2b;
  auto nccs = simplex_degree(cell_dim, cell_dim - 1);
  auto c2dim = mesh->get_array<I8>(cell_dim, "class_dim");
  auto s2dim = mesh->get_array<I8>(cell_dim - 1, "class_dim");
  auto ncands = cands2edges.size();
  auto cand_codes_w = Write<I8>(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto code = cand_codes[cand];
    auto e = cands2edges[cand];
    for (auto ec = e2ec[e]; ec < e2ec[e + 1]; ++ec) {
      auto c = ec2c[ec];
      auto ec_code = ec_codes[ec];
      auto cce = code_which_down(ec_code);
      auto rot = code_rotation(ec_code);
      auto c_dim = c2dim[c];
      for (Int eev_col = 0; eev_col < 2; ++eev_col) {
        if (!collapses(code, eev_col)) continue;
        auto eev_onto = 1 - eev_col;
        auto cev_onto = rot ^ eev_onto;
        auto ccv_onto = simplex_down_template(cell_dim, EDGE, cce, cev_onto);
        auto ccs_opp = simplex_opposite_template(cell_dim, VERT, ccv_onto);
        auto s_opp = cs2s[c * nccs + ccs_opp];
        if (s2dim[s_opp] != c_dim) {
          code = dont_collapse(code, eev_col);
        }
      }
    }
    cand_codes_w[cand] = code;
  };
  parallel_for(ncands, f, "check_collapse_exposure");
  return cand_codes_w;
}

Read<I8> check_collapse_exposure(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  for (Int cell_dim = EDGE + 1; cell_dim <= mesh->dim(); ++cell_dim) {
    cand_codes =
        check_collapse_exposure(mesh, cands2edges, cand_codes, cell_dim);
  }
  return mesh->sync_subset_array(
      EDGE, cand_codes, cands2edges, I8(DONT_COLLAPSE), 1);
}

}  // end namespace Omega_h
