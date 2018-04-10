#ifndef OMEGA_H_SIMPLIFY_HPP
#define OMEGA_H_SIMPLIFY_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

LOs tris_from_quads(LOs qv2v);
LOs tets_from_hexes(LOs hv2v);

void tris_from_quads_symmetric(Mesh* mesh);
void tets_from_hexes_symmetric(Mesh* mesh);

}  // end namespace Omega_h

#endif
