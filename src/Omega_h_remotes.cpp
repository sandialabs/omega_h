#include "Omega_h_remotes.hpp"

#include "Omega_h_map.hpp"

namespace Omega_h {

Remotes expand(Remotes a2c, LOs a2b) {
  return Remotes(expand(a2c.ranks, a2b, 1), expand(a2c.idxs, a2b, 1));
}

Remotes unmap(LOs a2b, Remotes b2c) {
  return Remotes(unmap(a2b, b2c.ranks, 1), unmap(a2b, b2c.idxs, 1));
}

Remotes identity_remotes(CommPtr comm, LO n) {
  return Remotes(Read<I32>(n, comm->rank()), LOs(n, 0, 1));
}

}  // end namespace Omega_h
