#ifndef SWAP3D_CHOICE_HPP
#define SWAP3D_CHOICE_HPP

#include "Omega_h_swap3d_loop.hpp"

namespace Omega_h {

namespace swap3d {

struct Choice {
  Int mesh;
  Real quality;
};

template <typename QualityMeasure, typename LengthMeasure>
OMEGA_H_DEVICE Choice choose(Loop loop, QualityMeasure const& quality_measure,
    LengthMeasure const& length_measure, Real max_length_allowed) {
  auto nmeshes = swap_mesh_counts[loop.size];
  auto nmesh_tris = swap_mesh_sizes[loop.size];
  auto uniq_tris2loop_verts = swap_triangles[loop.size];
  bool uniq_tris_cached[MAX_UNIQUE_TRIS] = {false};
  Real uniq_tri_quals[MAX_UNIQUE_TRIS] = {0};
  bool uniq_edgs_cached[MAX_UNIQUE_EDGES] = {false};
  Real uniq_edg_lens[MAX_UNIQUE_EDGES] = {0};
  Choice choice;
  choice.mesh = -1;
  choice.quality = 0.0;
  for (Int mesh = 0; mesh < nmeshes; ++mesh) {
    Real mesh_minqual = 1.0;
    auto mesh_tris2uniq_tris = &swap_meshes[loop.size][mesh * nmesh_tris];
    for (Int mesh_tri = 0; mesh_tri < nmesh_tris; ++mesh_tri) {
      auto uniq_tri = mesh_tris2uniq_tris[mesh_tri];
      if (!uniq_tris_cached[uniq_tri]) {
        auto tri_verts2loop_verts = uniq_tris2loop_verts[uniq_tri];
        /* the first three tet vertices are
           the same as the bottom triangle,
           curling into the tet. we fill these
           in from the triangle table for the current
           2D mesh being explored */
        Few<LO, 4> tet_verts2verts;
        for (Int tri_vert = 0; tri_vert < 3; ++tri_vert) {
          auto loop_vert = tri_verts2loop_verts[tri_vert];
          auto vert = loop.loop_verts2verts[loop_vert];
          tet_verts2verts[tri_vert] = vert;
        }
        /* each triangle will support two tets,
           one above and one below. this loop
           forms those tets, swapping vertices
           in between to maintain proper orientation.
           (mfr means Region of Face of Mesh) */
        Real tri_minqual = 1.0;
        for (Int tri_tet = 0; tri_tet < 2; ++tri_tet) {
          tet_verts2verts[3] = loop.eev2v[1 - tri_tet];
          auto tet_qual = quality_measure.measure(tet_verts2verts);
          tri_minqual = min2(tri_minqual, tet_qual);
          swap2(tet_verts2verts[1], tet_verts2verts[2]);
        }
        uniq_tris_cached[uniq_tri] = true;
        uniq_tri_quals[uniq_tri] = tri_minqual;
      }
      auto tri_minqual = uniq_tri_quals[uniq_tri];
      mesh_minqual = min2(mesh_minqual, tri_minqual);
      /* if we know this swap configuration will make
         negative tets, don't bother computing the rest of it. */
      if (mesh_minqual < 0.0) break;
    }
    if (mesh_minqual > choice.quality) {
      /* now that we have a configuration which is a candidate
         for being the new best, we will go ahead and check for
         edge length overshooting.
         the idea is to minimize the cost of this check.
         We'll use the same caching as we did for triangles above. */
      bool does_overshoot = false;
      for (Int mesh_edge = 0; mesh_edge < nedges[loop.size]; ++mesh_edge) {
        auto uniq_edge = edges2unique[loop.size][mesh][mesh_edge];
        if (!uniq_edgs_cached[uniq_edge]) {
          Few<LO, 2> edge_verts2verts;
          for (Int edge_vert = 0; edge_vert < 2; ++edge_vert) {
            auto loop_vert = unique_edges[loop.size][uniq_edge][edge_vert];
            auto vert = loop.loop_verts2verts[loop_vert];
            edge_verts2verts[edge_vert] = vert;
          }
          auto edge_length = length_measure.measure(edge_verts2verts);
          uniq_edg_lens[uniq_edge] = edge_length;
          uniq_edgs_cached[uniq_edge] = true;
        }
        auto edge_length = uniq_edg_lens[uniq_edge];
        if (edge_length > max_length_allowed) {
          does_overshoot = true;
          break;
        }
      }
      if (!does_overshoot) {
        choice.mesh = mesh;
        choice.quality = mesh_minqual;
      }
    }
  }
  return choice;
}

}  // end namespace swap3d

}  // end namespace Omega_h

#endif
