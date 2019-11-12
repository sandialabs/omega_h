#include "Omega_h_unmap_mesh.hpp"

#include "Omega_h_amr.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_profile.hpp"

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
  auto deg = element_degree(old_mesh->family(), ent_dim, ent_dim - 1);
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
  OMEGA_H_TIME_FUNCTION;
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

void unmap_parents(Mesh* old_mesh, Mesh* new_mesh,
    LOs new_ents2old_ents_a[], Few<LOs, 4> old_ents2new_ents) {
  for (Int ent_dim = 1; ent_dim <= new_mesh->dim(); ++ent_dim) {
    auto new_ents2old_ents = new_ents2old_ents_a[ent_dim];
    auto nnew_ents = new_ents2old_ents.size();
    Write<LO> new_parent_idx(nnew_ents, -1, "parent idx");
    Write<Byte> new_parent_code(nnew_ents, 0, "parent code");
    auto old_parents = old_mesh->ask_parents(ent_dim);
    auto functor = OMEGA_H_LAMBDA(LO new_ent) {
      auto old_ent = new_ents2old_ents[new_ent];
      auto old_parent_idx = old_parents.parent_idx[old_ent];
      auto old_parent_code = old_parents.codes[old_ent];
      if (old_parent_idx > -1) {
        auto old_parent_dim = amr::code_parent_dim(old_parent_code);
        new_parent_idx[new_ent] =
          old_ents2new_ents[old_parent_dim][old_parent_idx];
        if (new_parent_idx[new_ent] > -1)
          new_parent_code[new_ent] = old_parent_code;
      }
    };
    parallel_for(nnew_ents, std::move(functor));
    new_mesh->set_parents(ent_dim, Parents{new_parent_idx, new_parent_code});
  }
}

void unmap_leaves(Mesh* new_mesh) {
  for (Int ent_dim = 1; ent_dim <= new_mesh->dim(); ++ent_dim) {
    Write<Byte> leaf(new_mesh->nents(ent_dim));
    auto is_ent_leaf = new_mesh->ask_leaves(ent_dim);
    auto children = new_mesh->ask_children(ent_dim, ent_dim);
    auto functor = OMEGA_H_LAMBDA(LO ent) {
      leaf[ent] = is_ent_leaf[ent];
      auto nchild = children.a2ab[ent + 1] - children.a2ab[ent];
      if ((nchild == 0) && (!leaf[ent])) leaf[ent] = 1;
    };
    parallel_for(new_mesh->nents(ent_dim), std::move(functor));
    new_mesh->set_tag(ent_dim, "leaf", Omega_h::read(leaf));
  }
}

void unmap_mesh(Mesh* mesh, LOs new_ents2old_ents[]) {
  auto new_mesh = mesh->copy_meta();
  auto nnew_verts = (new_ents2old_ents[0].exists())
                        ? new_ents2old_ents[0].size()
                        : mesh->nverts();
  new_mesh.set_verts(nnew_verts);
  LOs old_lows2new_lows;
  Few<LOs, 4> old_ents2new_ents;
  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    if (ent_dim > VERT) {
      unmap_down(mesh, &new_mesh, ent_dim, new_ents2old_ents[ent_dim],
          old_lows2new_lows);
    }
    unmap_tags(mesh, &new_mesh, ent_dim, new_ents2old_ents[ent_dim]);
    old_ents2new_ents[ent_dim] =
        invert_injective_map(new_ents2old_ents[ent_dim], mesh->nents(ent_dim));
    unmap_owners(mesh, &new_mesh, ent_dim, new_ents2old_ents[ent_dim],
        old_ents2new_ents[ent_dim]);
    old_lows2new_lows = old_ents2new_ents[ent_dim];
  }
  if (mesh->has_any_parents()) {
    unmap_parents(
        mesh, &new_mesh, new_ents2old_ents, old_ents2new_ents);
    unmap_leaves(&new_mesh);
  }
  *mesh = new_mesh;
}

}  // end namespace Omega_h
