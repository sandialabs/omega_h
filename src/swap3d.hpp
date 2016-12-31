#ifndef SWAP3D_HPP
#define SWAP3D_HPP

#include "host_few.hpp"
#include "internal.hpp"

namespace Omega_h {

void swap3d_qualities(Mesh* mesh, AdaptOpts const& opts, LOs cands2edges,
    Reals* cand_quals, Read<I8>* cand_configs);

HostFew<LOs, 4> swap3d_keys_to_prods(Mesh* mesh, LOs keys2edges);

HostFew<LOs, 4> swap3d_topology(Mesh* mesh, LOs keys2edges,
    Read<I8> edge_configs, HostFew<LOs, 4> keys2prods);

bool swap_edges_3d(Mesh* mesh, AdaptOpts const& opts);

}  // end namespace Omega_h

#endif
