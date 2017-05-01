#ifndef OMEGA_H_BCAST_HPP
#define OMEGA_H_BCAST_HPP

#include <Omega_h_comm.hpp>

namespace Omega_h {

class Mesh;

void bcast_mesh(Mesh* mesh, CommPtr new_comm, bool is_source);

}  // end namespace Omega_h

#endif
