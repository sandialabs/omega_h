#include "Omega_h_motion.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_indset.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_modify.hpp"
#include "Omega_h_transfer.hpp"

#include <iostream>
#include "Omega_h_file.hpp"

namespace Omega_h {

static bool move_verts_ghosted(Mesh* mesh, AdaptOpts const& opts) {
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto comm = mesh->comm();
  auto verts_are_interior = mark_by_class_dim(mesh, VERT, mesh->dim());
  auto verts_are_cands = invert_marks(verts_are_interior);
  auto cands2verts = collect_marked(verts_are_cands);
  auto choices = get_motion_choices(mesh, opts, cands2verts);
  verts_are_cands =
      map_onto(choices.did_move, cands2verts, mesh->nverts(), I8(0), 1);
  if (get_sum(comm, verts_are_cands) == 0) return false;
  auto vert_quals =
      map_onto(choices.new_quals, cands2verts, mesh->nverts(), -1.0, 1);
  auto verts_are_keys = find_indset(mesh, VERT, vert_quals, verts_are_cands);
  mesh->add_tag(VERT, "key", 1, verts_are_keys);
  auto new_coords = deep_copy(mesh->coords());
  map_into(choices.new_coords, cands2verts, new_coords, mesh->dim());
  auto debug_old_coords = mesh->coords();
  for (LO i = 0; i < verts_are_keys.size(); ++i) {
    if (verts_are_keys[i]) {
      auto dist = norm(get_vector<3>(new_coords, i) - get_vector<3>(debug_old_coords, i));
      std::cout << "moving distance " << dist << '\n';
    }
  }
  mesh->add_tag(VERT, "motion_coords", mesh->dim(), Reals(new_coords));
  auto keys2verts = collect_marked(verts_are_keys);
  auto verts2cav_elems = mesh->ask_up(VERT, mesh->dim());
  set_owners_by_indset(mesh, VERT, keys2verts, verts2cav_elems);
  return true;
}

static void move_verts_elem_based(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto verts_are_keys = mesh->get_array<I8>(VERT, "key");
  mesh->remove_tag(VERT, "key");
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
    if (ent_dim == VERT) {
      new_mesh.set_verts(mesh->nverts());
    } else {
      new_mesh.set_ents(ent_dim, mesh->ask_down(ent_dim, ent_dim - 1));
    }
    new_mesh.set_owners(ent_dim, mesh->ask_owners(ent_dim));
    transfer_motion(
        mesh, opts.xfer_opts, &new_mesh, verts_are_keys, keys2verts, ent_dim);
    if (ent_dim == VERT) {
      auto new_coords = mesh->get_array<Real>(VERT, "motion_coords");
      mesh->remove_tag(VERT, "motion_coords");
      new_mesh.set_tag(VERT, "coordinates", new_coords);
    }
  }
  *mesh = new_mesh;
}

bool move_verts_to_conserve_size(Mesh* mesh, AdaptOpts const& opts) {
  if (!move_verts_ghosted(mesh, opts)) return false;
  mesh->set_parting(OMEGA_H_ELEM_BASED, false);
  move_verts_elem_based(mesh, opts);
  return true;
}

}  // end namespace Omega_h
