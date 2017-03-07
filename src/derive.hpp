#ifndef DERIVE_HPP
#define DERIVE_HPP

#include "Omega_h.hpp"

namespace Omega_h {

Reals derive_element_gradients(Mesh* mesh, Reals vert_values);
Reals derive_element_hessians(Mesh* mesh, Reals vert_gradients);
Reals recover_gradients(Mesh* mesh, Reals vert_values);
Reals recover_hessians_from_gradients(Mesh* mesh, Reals vert_gradients);
}  // namespace Omega_h

#endif
