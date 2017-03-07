#ifndef PROJECT_HPP
#define PROJECT_HPP

#include "Omega_h.hpp"

namespace Omega_h {

bool has_interior_verts(Mesh* mesh);
Reals project_by_fit(Mesh* mesh, Reals e_data);
Reals project_by_average(Mesh* mesh, Reals e_data);
}  // namespace Omega_h

#endif
