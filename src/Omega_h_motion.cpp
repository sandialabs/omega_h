#include "Omega_h_motion.hpp"

namespace Omega_h {

static bool motion_ghosted(Mesh* mesh, AdaptOpts const& opts) {
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto comm = mesh->comm();
  auto elems_are_cands =
      mark_sliver_layers(mesh, opts.min_quality_desired, opts.nsliver_layers);
  CHECK(get_max(comm, elems_are_cands) == 1);
  auto verts_are_cands = mark_down(mesh, mesh->dim(), VERT, elems_are_cands);
  auto cands2verts = collect_marked(edges_are_cands);
  auto choices = get_motion_choices(mesh, opts, cands2verts);
  verts_are_cands = choices.cands_did_move;
  if (sum(comm, verts_are_cands) == 0) return false;
  auto vert_quals = map_onto(choices.quals, cands2verts, mesh->nverts(), -1.0, 1);
  auto verts_are_keys = find_indset(mesh, VERT, vert_quals, verts_are_cands);
  auto new_coords_w = deep_copy(mesh->coords());
  map_into(choices.coords, cands2verts, new_coords_w, mesh->dim());
  mesh->add_tag(VERT, "motion_coords", mesh->dim(), OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DONT_OUTPUT, Reals(new_coords_w));
  auto new_metrics_w = deep_copy(mesh->coords());
  map_into(choices.coords, cands2verts, new_coords_w, mesh->dim());
  mesh->add_tag(VERT, "motion_coords", mesh->dim(), OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DONT_OUTPUT, Reals(new_coords_w));
  Graph edges2cav_elems;
  edges2cav_elems = mesh->ask_up(EDGE, mesh->dim());
  mesh->add_tag(EDGE, "config", 1, OMEGA_H_DONT_TRANSFER, OMEGA_H_DONT_OUTPUT,
      edge_configs);
  auto keys2edges = collect_marked(edges_are_keys);
  set_owners_by_indset(mesh, EDGE, keys2edges, edges2cav_elems);
  return true;
}

bool move_vertices_for_quality(Mesh* mesh, AdaptOpts const& opts) {
}

} // end namespace Omega_h
