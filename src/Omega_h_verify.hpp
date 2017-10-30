#ifndef OMEGA_H_VERIFY_HPP
#define OMEGA_H_VERIFY_HPP

namespace Omega_h {

class Mesh;

bool verify_down_verts(Mesh* mesh);
void verify_no_duplicates(Mesh* mesh);

}  // end namespace Omega_h

#endif
