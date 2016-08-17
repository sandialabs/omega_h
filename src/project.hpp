#ifndef PROJECT_HPP
#define PROJECT_HPP

#include "omega_h.hpp"

namespace Omega_h {

bool has_interior_verts(Mesh* mesh);
Reals project(Mesh* mesh, Reals e_data);
}

#endif
