#include "reorder.hpp"

#include "Omega_h_array_ops.hpp"
#include "hilbert.hpp"
#include "Omega_h_map.hpp"
#include "unmap_mesh.hpp"

namespace Omega_h {

/* Construct a graph from each vertex
   to the entities which use that vertex
   as their first vertex in the canonical
   downward ordering.
   The point of this is just to establish,
   for each entity, a single vertex that is
   "responsible for" that entity */
Graph find_entities_of_first_vertices(Mesh* mesh, Int ent_dim) {
  auto ev2v = mesh->ask_verts_of(ent_dim);
  auto e2fv = get_component(ev2v, ent_dim + 1, 0);
  auto fv2e = invert_map_by_sorting(e2fv, mesh->nverts());
  return fv2e;
}

LOs ent_order_from_vert_order(
    Mesh* mesh, Int ent_dim, LOs new_verts2old_verts) {
  CHECK(new_verts2old_verts.size() == mesh->nverts());
  auto old_verts2old_ents = find_entities_of_first_vertices(mesh, ent_dim);
  CHECK(old_verts2old_ents.a2ab.size() == mesh->nverts() + 1);
  CHECK(old_verts2old_ents.a2ab.last() == mesh->nents(ent_dim));
  CHECK(old_verts2old_ents.ab2b.size() == mesh->nents(ent_dim));
  auto new_verts2old_ents =
      unmap_graph(new_verts2old_verts, old_verts2old_ents);
  CHECK(new_verts2old_ents.a2ab.size() == mesh->nverts() + 1);
  CHECK(new_verts2old_ents.a2ab.last() == mesh->nents(ent_dim));
  CHECK(new_verts2old_ents.ab2b.size() == mesh->nents(ent_dim));
  auto new_ents2old_ents = new_verts2old_ents.ab2b;
  CHECK(new_ents2old_ents.size() == mesh->nents(ent_dim));
  return new_ents2old_ents;
}

void reorder_mesh(Mesh* old_mesh, Mesh* new_mesh, LOs new_verts2old_verts) {
  LOs new_ents2old_ents[4];
  new_ents2old_ents[VERT] = new_verts2old_verts;
  for (Int ent_dim = 1; ent_dim <= old_mesh->dim(); ++ent_dim) {
    new_ents2old_ents[ent_dim] =
        ent_order_from_vert_order(old_mesh, ent_dim, new_verts2old_verts);
  }
  unmap_mesh(old_mesh, new_mesh, new_ents2old_ents);
}

void reorder_mesh(Mesh* mesh, LOs new_verts2old_verts) {
  auto new_mesh = mesh->copy_meta();
  reorder_mesh(mesh, &new_mesh, new_verts2old_verts);
  *mesh = new_mesh;
}

void reorder_by_hilbert(Mesh* mesh) {
  auto coords = mesh->coords();
  LOs new_verts2old_verts = hilbert::sort_coords(coords, mesh->dim());
  reorder_mesh(mesh, new_verts2old_verts);
}

}  // end namespace Omega_h
