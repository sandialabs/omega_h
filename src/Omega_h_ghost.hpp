#ifndef OMEGA_H_GHOST_HPP
#define OMEGA_H_GHOST_HPP

#include <Omega_h_dist.hpp>

namespace Omega_h {

class Mesh;

/* a graph from local items to global items.
 * locals2edges is an offset map from local items to edges
 * outgoing from local items.
 * edges2remotes is a map from edges outgoing from local
 * items, to the remote destinations of said edges.
 */
struct RemoteGraph {
  LOs locals2edges;
  Remotes edges2remotes;
};

Dist get_local_elem_uses2own_verts(Mesh* mesh);
Remotes get_local_elem_uses2own_elems(Mesh* mesh);
void get_own_verts2own_elem_uses(
    Mesh* mesh, Remotes& serv_uses2own_elems, LOs& own_verts2serv_uses);
Remotes push_elem_uses(RemoteGraph own_verts2own_elems, Dist own_verts2verts);

void ghost_mesh(Mesh* mesh, Int nlayers, bool verbose);
void partition_by_verts(Mesh* mesh, bool verbose);
void partition_by_elems(Mesh* mesh, bool verbose);

}  // end namespace Omega_h

#endif
