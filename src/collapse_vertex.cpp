#include "coarsen.hpp"
#include "collapse.hpp"
#include "map.hpp"
#include "loop.hpp"

/* this file deals with the choice of which
 * edge collapse a vertex will participate in.
 */

namespace Omega_h {

void choose_vertex_collapses(Mesh* mesh, LOs cands2edges,
    Read<I8> cand_edge_codes, Reals cand_edge_quals, Read<I8>& verts_are_cands,
    Reals& vert_quals) {
  CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto edges2cands = invert_injective_map(cands2edges, mesh->nedges());
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ve_codes = v2e.codes;
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto verts_are_cands_w = Write<I8>(mesh->nverts());
  auto vert_quals_w = Write<Real>(mesh->nverts());
  auto f = LAMBDA(LO v) {
    bool vert_is_cand = false;
    Real best_qual = -1.0;
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      auto e = ve2e[ve];
      auto cand = edges2cands[e];
      if (cand == -1) continue;
      auto ve_code = ve_codes[ve];
      auto eev = code_which_down(ve_code);
      auto cand_code = cand_edge_codes[cand];
      if (!collapses(cand_code, eev)) continue;
      auto qual = cand_edge_quals[cand * 2 + eev];
      if (qual > best_qual) {
        vert_is_cand = true;
        best_qual = qual;
      }
    }
    verts_are_cands_w[v] = vert_is_cand;
    vert_quals_w[v] = best_qual;
  };
  parallel_for(mesh->nverts(), f);
  verts_are_cands = verts_are_cands_w;
  verts_are_cands = mesh->sync_array(VERT, verts_are_cands, 1);
  vert_quals = vert_quals_w;
  vert_quals = mesh->sync_array(VERT, vert_quals, 1);
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

void find_rails(Mesh* mesh, LOs keys2verts, Reals vert_quals,
    Read<I8> edge_cand_codes, Reals edge_cand_quals, LOs& rails2edges,
    Read<I8>& rail_col_dirs) {
  auto nkeys = keys2verts.size();
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ve_codes = v2e.codes;
  auto rails2edges_w = Write<LO>(nkeys, -1);
  auto rail_col_dirs_w = Write<I8>(nkeys, -1);
  auto f = LAMBDA(LO key) {
    auto v = keys2verts[key];
    auto vert_qual = vert_quals[v];
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      auto e = ve2e[ve];
      auto ve_code = ve_codes[ve];
      auto eev = code_which_down(ve_code);
      auto cand_code = edge_cand_codes[e];
      if (!collapses(cand_code, eev)) continue;
      auto edge_qual = edge_cand_quals[e * 2 + eev];
      if (edge_qual == vert_qual) {
        rails2edges_w[key] = e;
        rail_col_dirs_w[key] = static_cast<I8>(eev);
        return;
      }
    }
  };
  parallel_for(nkeys, f);
  rails2edges = rails2edges_w;
  rail_col_dirs = rail_col_dirs_w;
}

}
