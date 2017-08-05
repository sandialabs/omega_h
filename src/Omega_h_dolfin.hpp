#ifndef OMEGA_H_DOLFIN_HPP
#define OMEGA_H_DOLFIN_HPP

namespace dolfin {
class Mesh;
}

namespace Omega_h {

class Mesh;

void to_dolfin(dolfin::Mesh& mesh_dolfin, Mesh* mesh_osh);
void from_dolfin(Mesh* mesh_osh, dolfin::Mesh const& mesh_dolfin);

}

#endif // OMEGA_H_DOLFIN_HPP
