#ifndef INDSET_HPP
#define INDSET_HPP

#include "internal.hpp"

namespace Omega_h {

Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Reals quality, Read<I8> candidates);

}  // end namespace Omega_h

#endif
