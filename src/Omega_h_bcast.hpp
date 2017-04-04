#ifndef BCAST_HPP
#define BCAST_HPP

#include "Omega_h_internal.hpp"

namespace Omega_h {

void bcast_mesh(Mesh* mesh, CommPtr new_comm, bool is_source);

}  // end namespace Omega_h

#endif
