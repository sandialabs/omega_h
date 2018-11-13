#include "Omega_h_array_ops.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_hilbert.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_profile.hpp"
#include "Omega_h_sort.hpp"
#include "Omega_h_unmap_mesh.hpp"

namespace Omega_h {

/* Construct a graph from each vertex
   to the entities which use that vertex
   as their first vertex in the canonical
   downward ordering.
   The point of this is just to establish,
   for each entity, a single vertex that is
   "responsible for" that entity */
static Graph find_entities_of_first_vertices(Mesh* mesh, Int ent_dim) {
  auto ev2v = mesh->ask_verts_of(ent_dim);
  auto deg = element_degree(mesh->family(), ent_dim, VERT);
  auto e2fv = get_component(ev2v, deg, 0);
  auto fv2e = invert_map_by_sorting(e2fv, mesh->nverts());
  return fv2e;
}

static LOs ent_order_from_vert_order(
    Mesh* mesh, Int ent_dim, LOs new_verts2old_verts) {
  OMEGA_H_CHECK(new_verts2old_verts.size() == mesh->nverts());
  auto old_verts2old_ents = find_entities_of_first_vertices(mesh, ent_dim);
  OMEGA_H_CHECK(old_verts2old_ents.a2ab.size() == mesh->nverts() + 1);
  OMEGA_H_CHECK(old_verts2old_ents.a2ab.last() == mesh->nents(ent_dim));
  OMEGA_H_CHECK(old_verts2old_ents.ab2b.size() == mesh->nents(ent_dim));
  auto new_verts2old_ents =
      unmap_graph(new_verts2old_verts, old_verts2old_ents);
  OMEGA_H_CHECK(new_verts2old_ents.a2ab.size() == mesh->nverts() + 1);
  OMEGA_H_CHECK(new_verts2old_ents.a2ab.last() == mesh->nents(ent_dim));
  OMEGA_H_CHECK(new_verts2old_ents.ab2b.size() == mesh->nents(ent_dim));
  auto new_ents2old_ents = new_verts2old_ents.ab2b;
  OMEGA_H_CHECK(new_ents2old_ents.size() == mesh->nents(ent_dim));
  return new_ents2old_ents;
}

static void reorder_mesh_by_verts(Mesh* mesh, LOs new_verts2old_verts) {
  LOs new_ents2old_ents[4];
  new_ents2old_ents[VERT] = new_verts2old_verts;
  for (Int ent_dim = 1; ent_dim <= mesh->dim(); ++ent_dim) {
    new_ents2old_ents[ent_dim] =
        ent_order_from_vert_order(mesh, ent_dim, new_verts2old_verts);
  }
  unmap_mesh(mesh, new_ents2old_ents);
}

void reorder_by_hilbert(Mesh* mesh) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(mesh->comm()->size() == 1);
  auto coords = mesh->coords();
  LOs new_verts2old_verts = hilbert::sort_coords(coords, mesh->dim());
  reorder_mesh_by_verts(mesh, new_verts2old_verts);
  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    mesh->remove_tag(ent_dim, "global");
    mesh->add_tag(ent_dim, "global", 1, GOs(mesh->nents(ent_dim), 0, 1));
  }
}

void reorder_by_globals(Mesh* mesh) {
  LOs new_ents2old_ents[4];
  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    new_ents2old_ents[ent_dim] = sort_by_keys(mesh->globals(ent_dim));
  }
  unmap_mesh(mesh, new_ents2old_ents);
}

}  // end namespace Omega_h
