#include "Omega_h_motion.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_transfer.hpp"
#include "indset.hpp"
#include "modify.hpp"

#include <iostream>

namespace Omega_h {

static bool move_verts_ghosted(Mesh* mesh, AdaptOpts const& opts) {
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto comm = mesh->comm();
  auto elems_are_cands =
      mark_sliver_layers(mesh, opts.min_quality_desired, opts.nsliver_layers);
  CHECK(get_max(comm, elems_are_cands) == 1);
  auto verts_are_cands = mark_down(mesh, mesh->dim(), VERT, elems_are_cands);
  auto cands2verts = collect_marked(verts_are_cands);
  auto choices = get_motion_choices(mesh, opts, cands2verts);
  verts_are_cands =
      map_onto(choices.cands_did_move, cands2verts, mesh->nverts(), I8(0), 1);
  if (get_sum(comm, verts_are_cands) == 0) return false;
  auto vert_quals =
      map_onto(choices.quals, cands2verts, mesh->nverts(), -1.0, 1);
  auto verts_are_keys = find_indset(mesh, VERT, vert_quals, verts_are_cands);
  mesh->add_tag(VERT, "key", 1, verts_are_keys);
  auto ncomps = choices.new_sol.size() / mesh->nverts();
  mesh->add_tag(VERT, "motion_solution", ncomps, choices.new_sol);
  auto keys2verts = collect_marked(verts_are_keys);
  auto verts2cav_elems = mesh->ask_up(VERT, mesh->dim());
  set_owners_by_indset(mesh, VERT, keys2verts, verts2cav_elems);
  return true;
}

static void move_verts_elem_based(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto verts_are_keys = mesh->get_array<I8>(VERT, "key");
  mesh->remove_tag(VERT, "key");
  auto new_sol = mesh->get_array<Real>(VERT, "motion_solution");
  mesh->remove_tag(VERT, "motion_solution");
  auto keys2verts = collect_marked(verts_are_keys);
  if (opts.verbosity >= EACH_REBUILD) {
    auto nkeys = keys2verts.size();
    auto ntotal_keys = comm->allreduce(GO(nkeys), OMEGA_H_SUM);
    if (comm->rank() == 0) {
      std::cout << "moving " << ntotal_keys << " vertices\n";
    }
  }
  auto new_mesh = mesh->copy_meta();
  for (Int ent_dim = VERT; ent_dim <= mesh->dim(); ++ent_dim) {
    if (ent_dim == VERT)
      new_mesh.set_verts(mesh->nverts());
    else
      new_mesh.set_ents(ent_dim, mesh->ask_down(ent_dim, ent_dim - 1));
    new_mesh.set_owners(ent_dim, mesh->ask_owners(ent_dim));
    transfer_copy_motion(mesh, opts.xfer_opts, &new_mesh, ent_dim);
    if (ent_dim == VERT) {
      unpack_linearized_fields(
          mesh, opts.xfer_opts, &new_mesh, new_sol, verts_are_keys);
      if (mesh->has_tag(VERT, "warp")) {
        auto tb = mesh->get_tagbase(VERT, "warp");
        auto old_warp = mesh->get_array<Real>(VERT, "warp");
        auto old_coords = mesh->coords();
        auto new_coords = new_mesh.coords();
        auto motion = subtract_each(new_coords, old_coords);
        auto new_warp = subtract_each(old_warp, motion);
        new_mesh.add_tag(VERT, "warp", tb->ncomps(), new_warp);
      }
    } else if (ent_dim == EDGE) {
      auto edges_did_move = mark_up(&new_mesh, VERT, EDGE, verts_are_keys);
      auto new_edges2edges = collect_marked(edges_did_move);
      auto edges_didnt_move = invert_marks(edges_did_move);
      auto same_edges2edges = collect_marked(edges_didnt_move);
      transfer_length(
          mesh, &new_mesh, same_edges2edges, same_edges2edges, new_edges2edges);
    } else if (ent_dim == mesh->dim()) {
      auto elems_did_move =
          mark_up(&new_mesh, VERT, mesh->dim(), verts_are_keys);
      auto new_elems2elems = collect_marked(elems_did_move);
      auto elems_didnt_move = invert_marks(elems_did_move);
      auto same_elems2elems = collect_marked(elems_didnt_move);
      transfer_size(
          mesh, &new_mesh, same_elems2elems, same_elems2elems, new_elems2elems);
      transfer_quality(
          mesh, &new_mesh, same_elems2elems, same_elems2elems, new_elems2elems);
      auto verts2elems = mesh->ask_graph(VERT, mesh->dim());
      auto keys2elems = unmap_graph(keys2verts, verts2elems);
      transfer_pointwise(mesh, opts.xfer_opts, &new_mesh, VERT, keys2verts,
          keys2elems.a2ab, keys2elems.ab2b, same_elems2elems, same_elems2elems);
    }
  }
  *mesh = new_mesh;
}

bool move_verts_for_quality(Mesh* mesh, AdaptOpts const& opts) {
  if (!move_verts_ghosted(mesh, opts)) return false;
  mesh->set_parting(OMEGA_H_ELEM_BASED, false);
  move_verts_elem_based(mesh, opts);
  return true;
}

}  // end namespace Omega_h
