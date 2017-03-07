#include "refine.hpp"

#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "indset.hpp"
#include "Omega_h_map.hpp"
#include "modify.hpp"
#include "refine_qualities.hpp"
#include "refine_topology.hpp"
#include "transfer.hpp"

namespace Omega_h {

static bool refine_ghosted(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto edges_are_cands = mesh->get_array<I8>(EDGE, "candidate");
  mesh->remove_tag(EDGE, "candidate");
  auto cands2edges = collect_marked(edges_are_cands);
  auto cand_quals = refine_qualities(mesh, cands2edges);
  auto cands_are_good = each_geq_to(cand_quals, opts.min_quality_allowed);
  if (get_max(comm, cands_are_good) != 1) return false;
  auto nedges = mesh->nedges();
  auto edges_are_initial =
      map_onto(cands_are_good, cands2edges, nedges, I8(0), 1);
  auto edge_quals = map_onto(cand_quals, cands2edges, nedges, 0.0, 1);
  auto edges_are_keys = find_indset(mesh, EDGE, edge_quals, edges_are_initial);
  mesh->add_tag(EDGE, "key", 1, OMEGA_H_DONT_TRANSFER, OMEGA_H_DONT_OUTPUT,
      edges_are_keys);
  if (mesh->keeps_canonical_globals()) {
    mesh->add_tag(EDGE, "edge2rep_order", 1, OMEGA_H_DONT_TRANSFER,
        OMEGA_H_DONT_OUTPUT, get_edge2rep_order(mesh, edges_are_keys));
  }
  auto keys2edges = collect_marked(edges_are_keys);
  set_owners_by_indset(mesh, EDGE, keys2edges, mesh->ask_up(EDGE, mesh->dim()));
  return true;
}

static void refine_element_based(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto edges_are_keys = mesh->get_array<I8>(EDGE, "key");
  auto keys2edges = collect_marked(edges_are_keys);
  auto nkeys = keys2edges.size();
  auto ntotal_keys = comm->allreduce(GO(nkeys), OMEGA_H_SUM);
  if (opts.verbosity >= EACH_REBUILD && comm->rank() == 0) {
    std::cout << "refining " << ntotal_keys << " edges\n";
  }
  auto new_mesh = mesh->copy_meta();
  auto keys2midverts = LOs();
  auto old_verts2new_verts = LOs();
  auto old_lows2new_lows = LOs();
  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    auto keys2prods = LOs();
    auto prod_verts2verts = LOs();
    if (ent_dim == VERT) {
      keys2prods = LOs(nkeys + 1, 0, 1);
    } else {
      refine_products(mesh, ent_dim, keys2edges, keys2midverts,
          old_verts2new_verts, keys2prods, prod_verts2verts);
    }
    auto prods2new_ents = LOs();
    auto same_ents2old_ents = LOs();
    auto same_ents2new_ents = LOs();
    auto old_ents2new_ents = LOs();
    modify_ents(mesh, &new_mesh, ent_dim, EDGE, keys2edges, keys2prods,
        prod_verts2verts, old_lows2new_lows, &prods2new_ents,
        &same_ents2old_ents, &same_ents2new_ents, &old_ents2new_ents);
    if (ent_dim == VERT) {
      keys2midverts = prods2new_ents;
      old_verts2new_verts = old_ents2new_ents;
    }
    transfer_refine(mesh, &new_mesh, keys2edges, keys2midverts, ent_dim,
        keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents);
    old_lows2new_lows = old_ents2new_ents;
  }
  *mesh = new_mesh;
}

static bool refine(Mesh* mesh, AdaptOpts const& opts) {
  mesh->set_parting(OMEGA_H_GHOSTED);
  if (!refine_ghosted(mesh, opts)) return false;
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  refine_element_based(mesh, opts);
  return true;
}

bool refine_by_size(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto lengths = mesh->ask_lengths();
  auto edge_is_cand = each_gt(lengths, opts.max_length_desired);
  if (get_max(comm, edge_is_cand) != 1) return false;
  mesh->add_tag(EDGE, "candidate", 1, OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DONT_OUTPUT, edge_is_cand);
  return refine(mesh, opts);
}

}  // end namespace Omega_h
