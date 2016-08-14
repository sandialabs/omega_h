#ifndef DERIVE_HPP
#define DERIVE_HPP

#include "omega_h.hpp"

namespace osh {

Reals derive_element_gradients(Mesh* mesh, Reals vert_values);
Reals derive_element_hessians(Mesh* mesh, Reals vert_gradients);
Reals recover_gradients(Mesh* mesh, Reals vert_values);
Reals recover_hessians_from_gradients(Mesh* mesh, Reals vert_gradients);
Reals recover_hessians(Mesh* mesh, Reals vert_values);

}

#endif
