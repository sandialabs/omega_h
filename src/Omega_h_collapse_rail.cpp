#include "Omega_h_coarsen.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_collapse.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"

/* this file deals with the choice of which
 * edge collapse a vertex will participate in.
 */

namespace Omega_h {

void choose_rails(Mesh* mesh, LOs cands2edges, Read<I8> cand_edge_codes,
    Reals cand_edge_quals, Read<I8>* verts_are_cands, Reals* vert_quals,
    Read<GO>* vert_rails) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto edges2cands = invert_injective_map(cands2edges, mesh->nedges());
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ve_codes = v2e.codes;
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto lengths = mesh->ask_lengths();
  auto globals = mesh->globals(EDGE);
  auto verts_are_cands_w = Write<I8>(mesh->nverts());
  auto vert_quals_w = Write<Real>(mesh->nverts());
  auto vert_rails_w = Write<GO>(mesh->nverts());
  auto f = OMEGA_H_LAMBDA(LO v) {
    bool vert_is_cand = false;
    GO best_global = -1;
    Real best_length = -1;
    Real best_qual = -1;
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      auto e = ve2e[ve];
      auto cand = edges2cands[e];
      if (cand == -1) continue;
      auto ve_code = ve_codes[ve];
      auto eev = code_which_down(ve_code);
      auto cand_code = cand_edge_codes[cand];
      OMEGA_H_CHECK(cand_code != DONT_COLLAPSE);
      if (!collapses(cand_code, eev)) continue;
      auto global = globals[e];
      auto length = lengths[e];
      auto qual = cand_edge_quals[cand * 2 + eev];
      if ((best_global == -1) || (length < best_length) ||
          ((length == best_length) && (qual > best_qual))) {
        vert_is_cand = true;
        best_global = global;
        best_length = length;
        best_qual = qual;
      }
    }
    verts_are_cands_w[v] = vert_is_cand;
    if (vert_is_cand) OMEGA_H_CHECK(best_global != -1);
    vert_quals_w[v] = best_qual;
    vert_rails_w[v] = best_global;
  };
  parallel_for(mesh->nverts(), f, "choose_rails");
  *verts_are_cands = verts_are_cands_w;
  *verts_are_cands = mesh->sync_array(VERT, *verts_are_cands, 1);
  *vert_quals = vert_quals_w;
  *vert_quals = mesh->sync_array(VERT, *vert_quals, 1);
  *vert_rails = vert_rails_w;
  *vert_rails = mesh->sync_array(VERT, *vert_rails, 1);
}

/* this function is in some sense
   the inverse of choose_vertex_collapses(),
   and shares much the same logic.
   after the independent set of key vertices
   is chosen, we go back and figure out which
   of the adjacent edges it was collapsing across.
   we do this by matching quality, so if two
   edges were equally good, the first upward
   adjacent will be chosen */

void find_rails(Mesh* mesh, LOs keys2verts, Read<GO> verts2rail,
    LOs* rails2edges, Read<I8>* rail_col_dirs) {
  auto nkeys = keys2verts.size();
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ve_codes = v2e.codes;
  auto rails2edges_w = Write<LO>(nkeys, -1);
  auto rail_col_dirs_w = Write<I8>(nkeys, -1);
  auto edge_globals = mesh->globals(EDGE);
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto v = keys2verts[key];
    auto rail_global = verts2rail[v];
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      auto e = ve2e[ve];
      auto edge_global = edge_globals[e];
      if (edge_global == rail_global) {
        auto ve_code = ve_codes[ve];
        auto eev = code_which_down(ve_code);
        rails2edges_w[key] = e;
        rail_col_dirs_w[key] = static_cast<I8>(eev);
        return;
      }
    }
    OMEGA_H_NORETURN();
  };
  parallel_for(nkeys, f, "find_rails");
  *rails2edges = rails2edges_w;
  *rail_col_dirs = rail_col_dirs_w;
}
}  // namespace Omega_h
