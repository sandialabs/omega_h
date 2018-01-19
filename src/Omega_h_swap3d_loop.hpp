#ifndef SWAP3D_LOOP_HPP
#define SWAP3D_LOOP_HPP

#include "Omega_h_swap3d_tables.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_few.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

namespace swap3d {

/* by definition, the loop vertices curl
   around the edge by the right-hand rule,
   i.e. counterclockwise when looking from
   the second edge vertex to the first. */
struct Loop {
  Int size;
  Few<LO, 2> eev2v;
  Few<LO, MAX_EDGE_SWAP> loop_verts2verts;
};

OMEGA_H_DEVICE Loop find_loop(LOs const& edges2edge_tets,
    LOs const& edge_tets2tets, Read<I8> const& edge_tet_codes,
    LOs const& edge_verts2verts, LOs const& tet_verts2verts, LO edge) {
  Loop loop;
  auto begin_use = edges2edge_tets[edge];
  auto end_use = edges2edge_tets[edge + 1];
  loop.size = end_use - begin_use;
  if (loop.size > MAX_EDGE_SWAP) return loop;
  OMEGA_H_CHECK(loop.size >= 3);
  for (Int eev = 0; eev < 2; ++eev) {
    loop.eev2v[eev] = edge_verts2verts[edge * 2 + eev];
  }
  /* collect the endpoints of the loop edges.
     each pair of endpoints is chosen to be pointing
     in the direction of curl. */
  Few<LO, 2> tmp_edges[MAX_EDGE_SWAP];
  for (Int i = 0; i < MAX_EDGE_SWAP; ++i) {
    tmp_edges[i][0] = tmp_edges[i][1] = -1;
  }
  for (Int loop_edge = 0; loop_edge < loop.size; ++loop_edge) {
    auto edge_tet = begin_use + loop_edge;
    auto tet = edge_tets2tets[edge_tet];
    auto code = edge_tet_codes[edge_tet];
    auto rre = code_which_down(code);
    auto rot = code_rotation(code);
    auto rre_opp = simplex_opposite_template(REGION, EDGE, rre);
    for (Int eev = 0; eev < 2; ++eev) {
      /* we rely on the fact that tet-edge-vertices are
         defined such that the opposite edge curls around
         the input edge. */
      auto rev = rot ^ eev;
      auto rrv = simplex_down_template(REGION, EDGE, rre_opp, rev);
      auto v = tet_verts2verts[tet * 4 + rrv];
      tmp_edges[edge_tet - begin_use][eev] = v;
    }
  }
  /* Upward adjacency from edges to tets is not ordered
   * to curl around; it is ordered by global numbers.
   * The following code uses insertion sort to
   * order the edges around the loop by matching their
   * endpoints.
   * Remember, there are at most 7 edges to sort. */
  for (Int i = 0; i < loop.size - 1; ++i) {
    Int j;
    for (j = i + 1; j < loop.size; ++j) {
      if (tmp_edges[j][0] == tmp_edges[i][1]) break;
    }
    OMEGA_H_CHECK(j < loop.size);
    swap2(tmp_edges[i + 1], tmp_edges[j]);
  }
  for (Int lv = 0; lv < loop.size; ++lv) {
    loop.loop_verts2verts[lv] = tmp_edges[lv][0];
  }
  return loop;
}

}  // end namespace swap3d

}  // end namespace Omega_h

#endif
