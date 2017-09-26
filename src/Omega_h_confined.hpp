#ifndef OMEGA_H_CONFINED_HPP
#define OMEGA_H_CONFINED_HPP

#include "Omega_h_array.hpp"

namespace Omega_h {

class Mesh;

Bytes find_bridge_edges(Mesh* mesh);

Reals get_pad_dists(Mesh* mesh, Int pad_dim, Read<I8> edges_are_bridges);
Reals get_pinched_angles(Mesh* mesh, Int pinched_dim);

}  // end namespace Omega_h

#endif
