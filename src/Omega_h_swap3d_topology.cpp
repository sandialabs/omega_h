#include "Omega_h_swap3d.hpp"

#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_swap3d_loop.hpp"
#include "Omega_h_swap3d_tables.hpp"

namespace Omega_h {

HostFew<LOs, 4> swap3d_keys_to_prods(Mesh* mesh, LOs keys2edges) {
  auto edges2tets = mesh->ask_up(EDGE, REGION);
  auto edges2edge_tets = edges2tets.a2ab;
  auto edges2ntets = get_degrees(edges2edge_tets);
  auto nkeys = keys2edges.size();
  HostFew<Write<LO>, 4> keys2nprods_w;
  for (Int prod_dim = EDGE; prod_dim <= REGION; ++prod_dim) {
    keys2nprods_w[prod_dim] = Write<LO>(nkeys);
  }
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto edge = keys2edges[key];
    auto loop_size = edges2ntets[edge];
    auto nplane_tris = swap3d::swap_mesh_sizes[loop_size];
    auto nplane_edges = swap3d::nedges[loop_size];
    auto nprod_edges = nplane_edges;
    auto nprod_tris = nplane_tris + 2 * nplane_edges;
    auto nprod_tets = 2 * nplane_tris;
    keys2nprods_w[EDGE][key] = nprod_edges;
    keys2nprods_w[FACE][key] = nprod_tris;
    keys2nprods_w[REGION][key] = nprod_tets;
  };
  parallel_for(nkeys, f, "swap3d_keys_to_prods");
  HostFew<LOs, 4> keys2prods;
  for (Int prod_dim = EDGE; prod_dim <= REGION; ++prod_dim) {
    keys2prods[prod_dim] = offset_scan(LOs(keys2nprods_w[prod_dim]));
  }
  return keys2prods;
}

HostFew<LOs, 4> swap3d_topology(Mesh* mesh, LOs keys2edges,
    Read<I8> edge_configs, HostFew<LOs, 4> keys2prods) {
  auto edges2tets = mesh->ask_up(EDGE, REGION);
  auto edges2edge_tets = edges2tets.a2ab;
  auto edge_tets2tets = edges2tets.ab2b;
  auto edge_tet_codes = edges2tets.codes;
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto tet_verts2verts = mesh->ask_verts_of(REGION);
  HostFew<Write<LO>, 4> prod_verts2verts_w;
  for (Int prod_dim = EDGE; prod_dim <= REGION; ++prod_dim) {
    prod_verts2verts_w[prod_dim] =
        Write<LO>(keys2prods[prod_dim].last() * Int(prod_dim + 1));
  }
  auto nkeys = keys2edges.size();
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto edge = keys2edges[key];
    auto config = edge_configs[edge];
    auto loop = swap3d::find_loop(edges2edge_tets, edge_tets2tets,
        edge_tet_codes, edge_verts2verts, tet_verts2verts, edge);
    auto nplane_tris = swap3d::swap_mesh_sizes[loop.size];
    auto nplane_edges = swap3d::nedges[loop.size];
    for (Int plane_edge = 0; plane_edge < nplane_edges; ++plane_edge) {
      auto unique_edge = swap3d::edges2unique[loop.size][config][plane_edge];
      Few<LO, 2> plane_edge_verts;
      for (Int pev = 0; pev < 2; ++pev) {
        auto loop_vert = swap3d::unique_edges[loop.size][unique_edge][pev];
        auto vert = loop.loop_verts2verts[loop_vert];
        plane_edge_verts[pev] = vert;
      }
      auto prod_edge = keys2prods[EDGE][key] + plane_edge;
      for (Int pev = 0; pev < 2; ++pev) {
        prod_verts2verts_w[EDGE][prod_edge * 2 + pev] = plane_edge_verts[pev];
      }
      for (Int eev = 0; eev < 2; ++eev) {
        Few<LO, 3> new_tri_verts;
        for (Int nfv = 0; nfv < 2; ++nfv) {
          new_tri_verts[nfv] = plane_edge_verts[nfv];
        }
        new_tri_verts[2] = loop.eev2v[eev];
        auto prod_tri = keys2prods[FACE][key] + 2 * plane_edge + eev;
        for (Int nfv = 0; nfv < 3; ++nfv) {
          prod_verts2verts_w[FACE][prod_tri * 3 + nfv] = new_tri_verts[nfv];
        }
      }
    }
    for (Int plane_tri = 0; plane_tri < nplane_tris; ++plane_tri) {
      auto uniq_tri =
          swap3d::swap_meshes[loop.size][config * nplane_tris + plane_tri];
      Few<LO, 3> plane_tri_verts;
      for (Int pfv = 0; pfv < 3; ++pfv) {
        auto loop_vert = swap3d::swap_triangles[loop.size][uniq_tri][pfv];
        auto vert = loop.loop_verts2verts[loop_vert];
        plane_tri_verts[pfv] = vert;
      }
      auto prod_tri = keys2prods[FACE][key] + 2 * nplane_edges + plane_tri;
      for (Int pfv = 0; pfv < 3; ++pfv) {
        prod_verts2verts_w[FACE][prod_tri * 3 + pfv] = plane_tri_verts[pfv];
      }
      for (Int eev = 0; eev < 2; ++eev) {
        auto prod_tet = keys2prods[REGION][key] + 2 * plane_tri + eev;
        Few<LO, 4> new_tet_verts;
        for (Int pfv = 0; pfv < 3; ++pfv) {
          new_tet_verts[pfv] = plane_tri_verts[pfv];
        }
        if (eev == 0) swap2(new_tet_verts[1], new_tet_verts[2]);
        new_tet_verts[3] = loop.eev2v[eev];
        for (Int nrv = 0; nrv < 4; ++nrv) {
          prod_verts2verts_w[REGION][prod_tet * 4 + nrv] = new_tet_verts[nrv];
        }
      }
    }
  };
  parallel_for(nkeys, f, "swap3d_topology");
  HostFew<LOs, 4> prod_verts2verts;
  for (Int prod_dim = EDGE; prod_dim <= REGION; ++prod_dim) {
    prod_verts2verts[prod_dim] = prod_verts2verts_w[prod_dim];
  }
  return prod_verts2verts;
}

}  // end namespace Omega_h
