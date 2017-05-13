#include "Omega_h_unmap_mesh.hpp"

#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

void unmap_tags(
    Mesh* old_mesh, Mesh* new_mesh, Int ent_dim, LOs new_ents2old_ents) {
  for (Int i = 0; i < old_mesh->ntags(ent_dim); ++i) {
    auto tag = old_mesh->get_tag(ent_dim, i);
    if (is<I8>(tag)) {
      new_mesh->add_tag<I8>(ent_dim, tag->name(), tag->ncomps(),
          unmap(new_ents2old_ents, as<I8>(tag)->array(), tag->ncomps()));
    } else if (is<I32>(tag)) {
      new_mesh->add_tag<I32>(ent_dim, tag->name(), tag->ncomps(),
          unmap(new_ents2old_ents, as<I32>(tag)->array(), tag->ncomps()));
    } else if (is<I64>(tag)) {
      new_mesh->add_tag<I64>(ent_dim, tag->name(), tag->ncomps(),
          unmap(new_ents2old_ents, as<I64>(tag)->array(), tag->ncomps()));
    } else if (is<Real>(tag)) {
      new_mesh->add_tag<Real>(ent_dim, tag->name(), tag->ncomps(),
          unmap(new_ents2old_ents, as<Real>(tag)->array(), tag->ncomps()));
    }
  }
}

void unmap_down(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs new_ents2old_ents, LOs old_lows2new_lows) {
  auto deg = simplex_degrees[ent_dim][ent_dim - 1];
  auto old_ents2old_lows = old_mesh->ask_down(ent_dim, ent_dim - 1);
  auto oel2ol = old_ents2old_lows.ab2b;
  auto oe2l_codes = old_ents2old_lows.codes;
  auto nel2ol = unmap(new_ents2old_ents, oel2ol, deg);
  auto nel2nl = compound_maps(nel2ol, old_lows2new_lows);
  auto new_ents2new_lows = Adj(nel2nl);
  if (oe2l_codes.exists()) {
    auto ne2l_codes = unmap(new_ents2old_ents, oe2l_codes, deg);
    new_ents2new_lows.codes = ne2l_codes;
  }
  new_mesh->set_ents(ent_dim, new_ents2new_lows);
}

Remotes unmap_owners(
    Mesh* old_mesh, Int ent_dim, LOs new_ents2old_ents, LOs old_ents2new_ents) {
  auto old_copies2old_owners = old_mesh->ask_dist(ent_dim);
  auto old_owners2old_copies = old_copies2old_owners.invert();
  auto old_copies2new_owners = old_owners2old_copies.exch(old_ents2new_ents, 1);
  auto new_ents2new_owners = unmap(new_ents2old_ents, old_copies2new_owners, 1);
  auto old_own_ranks = old_mesh->ask_owners(ent_dim).ranks;
  auto new_own_ranks = unmap(new_ents2old_ents, old_own_ranks, 1);
  return Remotes(new_own_ranks, new_ents2new_owners);
}

void unmap_owners(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs new_ents2old_ents, LOs old_ents2new_ents) {
  if (old_mesh->comm()->size() == 1) return;
  auto owners =
      unmap_owners(old_mesh, ent_dim, new_ents2old_ents, old_ents2new_ents);
  new_mesh->set_owners(ent_dim, owners);
}

void unmap_mesh(Mesh* mesh, LOs new_ents2old_ents[]) {
  auto new_mesh = mesh->copy_meta();
  new_mesh.set_verts(mesh->nverts());
  LOs old_lows2new_lows;
  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    if (ent_dim > VERT) {
      unmap_down(mesh, &new_mesh, ent_dim, new_ents2old_ents[ent_dim],
          old_lows2new_lows);
    }
    unmap_tags(mesh, &new_mesh, ent_dim, new_ents2old_ents[ent_dim]);
    auto old_ents2new_ents = invert_injective_map(
        new_ents2old_ents[ent_dim], mesh->nents(ent_dim));
    unmap_owners(mesh, &new_mesh, ent_dim, new_ents2old_ents[ent_dim],
        old_ents2new_ents);
    old_lows2new_lows = old_ents2new_ents;
  }
  *mesh = new_mesh;
}

}  // end namespace Omega_h
