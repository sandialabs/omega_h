#ifndef GHOST_HPP
#define GHOST_HPP

#include "internal.hpp"

namespace Omega_h {

Dist get_local_elem_uses2own_verts(Mesh* mesh);
Remotes get_local_elem_uses2own_elems(Mesh* mesh);
void get_own_verts2own_elem_uses(
    Mesh* mesh, Remotes& serv_uses2own_elems, LOs& own_verts2serv_uses);
Remotes push_elem_uses(
    Remotes serv_uses2own_elems, LOs own_verts2serv_uses, Dist own_verts2verts);

void ghost_mesh(Mesh* mesh, bool verbose);
void partition_by_verts(Mesh* mesh, bool verbose);
void partition_by_elems(Mesh* mesh, bool verbose);

}  // end namespace Omega_h

#endif
