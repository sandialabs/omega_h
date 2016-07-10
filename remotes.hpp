#ifndef REMOTES_HPP
#define REMOTES_HPP

#include "internal.hpp"

namespace osh {

Remotes expand(Remotes a2c, LOs a2b);
Remotes unmap(LOs a2b, Remotes b2c);
Remotes identity_remotes(CommPtr comm, LO n);

} //end namespace osh

#endif
