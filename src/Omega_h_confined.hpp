#ifndef OMEGA_H_CONFINED_HPP
#define OMEGA_H_CONFINED_HPP

#include "Omega_h.hpp"

namespace Omega_h {
Read<I8> find_bridge_edges(Mesh* mesh);
Read<I8> find_angle_triangles(Mesh* mesh);
Read<I8> find_angle_elems(Mesh* mesh);
}

#endif
