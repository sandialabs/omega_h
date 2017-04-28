#ifndef OMEGA_H_REFINE_QUALITIES_HPP
#define OMEGA_H_REFINE_QUALITIES_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

Reals refine_qualities(Mesh* mesh, LOs candidates);

}  // end namespace Omega_h

#endif
