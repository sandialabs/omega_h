#ifndef OMEGA_H_SIMPLIFY_HPP
#define OMEGA_H_SIMPLIFY_HPP

namespace Omega_h {

class Mesh;

void quads2tris(Mesh* mesh);
void hexes2tets(Mesh* mesh);

}  // end namespace Omega_h

#endif
