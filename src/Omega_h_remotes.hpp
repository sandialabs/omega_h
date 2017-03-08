#ifndef OMEGA_H_REMOTES_HPP
#define OMEGA_H_REMOTES_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_comm.hpp>

namespace Omega_h {

struct Remotes {
  Remotes() {}
  Remotes(Read<I32> ranks_, LOs idxs_) : ranks(ranks_), idxs(idxs_) {}
  Read<I32> ranks;
  LOs idxs;
};

Remotes expand(Remotes a2c, LOs a2b);
Remotes unmap(LOs a2b, Remotes b2c);
Remotes identity_remotes(CommPtr comm, LO n);

Remotes owners_from_globals(
    CommPtr comm, Read<GO> globals, Read<I32> own_ranks);

}  // end namespace Omega_h

#endif
