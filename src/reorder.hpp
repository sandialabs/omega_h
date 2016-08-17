#ifndef REORDER_HPP
#define REORDER_HPP

#include "internal.hpp"

namespace Omega_h {

Graph find_entities_of_first_vertices(Mesh* mesh, Int ent_dim);
LOs ent_order_from_vert_order(Mesh* mesh, Int ent_dim, LOs new_vert2old_vert);
void reorder_mesh(Mesh* old_mesh, Mesh* new_mesh, LOs new_verts2old_verts);
void reorder_mesh(Mesh* mesh, LOs new_verts2old_verts);
void reorder_by_hilbert(Mesh* mesh);

}  // end namespace Omega_h

#endif
