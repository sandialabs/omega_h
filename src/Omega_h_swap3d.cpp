#include "Omega_h_swap3d.hpp"

#include <iostream>

#include "Omega_h_indset.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_modify.hpp"
#include "Omega_h_swap.hpp"
#include "Omega_h_transfer.hpp"

namespace Omega_h {

static bool swap3d_ghosted(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto edges_are_cands = mesh->get_array<I8>(EDGE, "candidate");
  mesh->remove_tag(EDGE, "candidate");
  auto cands2edges = collect_marked(edges_are_cands);
  auto cand_quals = Reals();
  auto cand_configs = Read<I8>();
  swap3d_qualities(mesh, opts, cands2edges, &cand_quals, &cand_configs);
  auto edge_configs =
      map_onto(cand_configs, cands2edges, mesh->nedges(), I8(-1), 1);
  auto keep_cands = filter_swap_improve(mesh, cands2edges, cand_quals);
  filter_swap(keep_cands, &cands2edges, &cand_quals);
  if (comm->reduce_and(cands2edges.size() == 0)) return false;
  edges_are_cands = mark_image(cands2edges, mesh->nedges());
  auto edge_quals = map_onto(cand_quals, cands2edges, mesh->nedges(), -1.0, 1);
  auto edges_are_keys = find_indset(mesh, EDGE, edge_quals, edges_are_cands);
  Graph edges2cav_elems;
  edges2cav_elems = mesh->ask_up(EDGE, mesh->dim());
  mesh->add_tag(EDGE, "key", 1, edges_are_keys);
  mesh->add_tag(EDGE, "config", 1, edge_configs);
  auto keys2edges = collect_marked(edges_are_keys);
  set_owners_by_indset(mesh, EDGE, keys2edges, edges2cav_elems);
  return true;
}

static void swap3d_element_based(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto edges_are_keys = mesh->get_array<I8>(EDGE, "key");
  mesh->remove_tag(EDGE, "key");
  auto edges_configs = mesh->get_array<I8>(EDGE, "config");
  mesh->remove_tag(EDGE, "config");
  auto keys2edges = collect_marked(edges_are_keys);
  if (opts.verbosity >= EACH_REBUILD) {
    auto nkeys = keys2edges.size();
    auto ntotal_keys = comm->allreduce(GO(nkeys), OMEGA_H_SUM);
    if (comm->rank() == 0) {
      std::cout << "swapping " << ntotal_keys << " 3D edges\n";
    }
  }
  auto new_mesh = mesh->copy_meta();
  new_mesh.set_verts(mesh->nverts());
  new_mesh.set_owners(VERT, mesh->ask_owners(VERT));
  transfer_copy(mesh, opts.xfer_opts, &new_mesh, VERT);
  auto keys2prods = swap3d_keys_to_prods(mesh, keys2edges);
  auto prod_verts2verts =
      swap3d_topology(mesh, keys2edges, edges_configs, keys2prods);
  auto old_lows2new_lows = LOs(mesh->nverts(), 0, 1);
  for (Int ent_dim = EDGE; ent_dim <= mesh->dim(); ++ent_dim) {
    auto prods2new_ents = LOs();
    auto same_ents2old_ents = LOs();
    auto same_ents2new_ents = LOs();
    auto old_ents2new_ents = LOs();
    modify_ents_adapt(mesh, &new_mesh, ent_dim, EDGE, keys2edges,
        keys2prods[ent_dim], prod_verts2verts[ent_dim], old_lows2new_lows,
        &prods2new_ents, &same_ents2old_ents, &same_ents2new_ents,
        &old_ents2new_ents);
    transfer_swap(mesh, opts.xfer_opts, &new_mesh, ent_dim, keys2edges,
        keys2prods[ent_dim], prods2new_ents, same_ents2old_ents,
        same_ents2new_ents);
    old_lows2new_lows = old_ents2new_ents;
  }
  *mesh = new_mesh;
}

bool swap_edges_3d(Mesh* mesh, AdaptOpts const& opts) {
  if (!swap_part1(mesh, opts)) return false;
  if (!swap3d_ghosted(mesh, opts)) return false;
  mesh->set_parting(OMEGA_H_ELEM_BASED, false);
  swap3d_element_based(mesh, opts);
  return true;
}

}  // end namespace Omega_h
