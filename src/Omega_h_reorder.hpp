#ifndef OMEGA_H_REORDER_HPP
#define OMEGA_H_REORDER_HPP

#include <Omega_h_graph.hpp>

namespace Omega_h {

class Mesh;

void reorder_by_hilbert(Mesh* mesh);

}  // end namespace Omega_h

#endif
