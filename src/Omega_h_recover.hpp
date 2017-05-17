#ifndef OMEGA_H_RECOVER_HPP
#define OMEGA_H_RECOVER_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

bool has_interior_verts(Mesh* mesh);
Reals project_by_fit(Mesh* mesh, Reals e_data);
Reals project_by_average(Mesh* mesh, Reals e_data);

Reals derive_element_gradients(Mesh* mesh, Reals vert_values);
Reals derive_element_hessians(Mesh* mesh, Reals vert_gradients);

Reals recover_gradients(Mesh* mesh, Reals vert_values);
Reals recover_hessians_from_gradients(Mesh* mesh, Reals vert_gradients);
Reals recover_hessians(Mesh* mesh, Reals vert_values);

}  // namespace Omega_h

#endif
