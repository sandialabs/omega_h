#ifndef OMEGA_H_UNMAP_MESH_HPP
#define OMEGA_H_UNMAP_MESH_HPP

#include <Omega_h_remotes.hpp>

namespace Omega_h {

class Mesh;

void unmap_tags(
    Mesh* old_mesh, Mesh* new_mesh, Int ent_dim, LOs new_ents2old_ents);

void unmap_down(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs new_ents2old_ents, LOs old_lows2new_lows);

/* Given an injective map (new_ents2old_ents) and its inverse
   (old_ents2new_ents), get the owners of (ent_dim) in the old mesh and unmap to
   what they will be in the new mesh. For each *owner*, new_index =
   old_ents2new_ents[old_index]. These new indices are then communicated to the
   copies and returned in the order indicated by (new_ents2old_ents) */
Remotes unmap_owners(
    Mesh* old_mesh, Int ent_dim, LOs new_ents2old_ents, LOs old_ents2new_ents);

/* Sets the owners of the new mesh based on the above */
void unmap_owners(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs new_ents2old_ents, LOs old_ents2new_ents);

void unmap_mesh(Mesh* mesh, LOs new_ents2old_ents[]);

}  // end namespace Omega_h

#endif
