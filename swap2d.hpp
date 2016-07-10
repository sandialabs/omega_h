#ifndef SWAP2D_HPP
#define SWAP2D_HPP

#include "internal.hpp"
#include "few.hpp"

namespace osh {

Reals swap2d_qualities(Mesh* mesh, LOs cands2edges);

void swap2d_topology(Mesh* mesh, LOs keys2edges,
                     HostFew<LOs, 3>* keys2prods_out,
                     HostFew<LOs, 3>* prod_verts2verts_out);

bool swap2d(Mesh* mesh, Real qual_ceil, Int nlayers, bool verbose);

} //end namespace osh

#endif
