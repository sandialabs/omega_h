#include "Omega_h_motion.hpp"

namespace Omega_h {

static bool move_verts_ghosted(Mesh* mesh, AdaptOpts const& opts) {
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto comm = mesh->comm();
  auto elems_are_cands =
      mark_sliver_layers(mesh, opts.min_quality_desired, opts.nsliver_layers);
  CHECK(get_max(comm, elems_are_cands) == 1);
  auto verts_are_cands = mark_down(mesh, mesh->dim(), VERT, elems_are_cands);
  auto cands2verts = collect_marked(edges_are_cands);
  auto choices = get_motion_choices(mesh, opts, cands2verts);
  verts_are_cands = map_onto(choices.cands_did_move, cands2verts, mesh->nverts(),
      I8(0), 1);
  if (sum(comm, verts_are_cands) == 0) return false;
  auto vert_quals = map_onto(choices.quals, cands2verts, mesh->nverts(), -1.0, 1);
  auto verts_are_keys = find_indset(mesh, VERT, vert_quals, verts_are_cands);
  mesh->add_tag(VERT, "key", 1, OMEGA_H_DONT_TRANSFER, OMEGA_H_DONT_OUTPUT,
      verts_are_keys);
  auto ncomps = choices.new_sol.size() / mesh->nverts();
  mesh->add_tag(VERT, "motion_solution", ncomps, OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DONT_OUTPUT, choices.new_sol);
  auto keys2verts = collect_marked(verts_are_keys);
  auto verts2cav_elems = mesh->ask_up(VERT, mesh->dim());
  set_owners_by_indset(mesh, VERT, keys2verts, verts2cav_elems);
  return true;
}

static void move_verts_element_based(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto verts_are_keys = mesh->get_array<I8>(VERT, "key");
  mesh->remove_tag(VERT, "key");
  auto new_sol = mesh->get_array<I8>(VERT, "motion_solution");
  mesh->remove_tag(VERT, "motion_solution");
  auto keys2verts = collect_marked(verts_are_keys);
  if (opts.verbosity >= EACH_REBUILD) {
    auto nkeys = keys2verts.size();
    auto ntotal_keys = comm->allreduce(GO(nkeys), OMEGA_H_SUM);
    if (comm->rank() == 0) {
      std::cout << "moving " << ntotal_keys << " vertices\n";
    }
  }
}

bool move_verts_for_quality(Mesh* mesh, AdaptOpts const& opts) {
  if (!move_verts_ghosted(mesh, opts)) return false;
  mesh->set_parting(OMEGA_H_ELEM_BASED, false);
  move_verts_elem_based(mesh, opts);
}

} // end namespace Omega_h
