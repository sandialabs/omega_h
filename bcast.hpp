#ifndef BCAST_HPP
#define BCAST_HPP

#include "internal.hpp"

namespace osh {

void bcast_mesh(Mesh* mesh, CommPtr new_comm, bool is_source);

} //end namespace osh

#endif
