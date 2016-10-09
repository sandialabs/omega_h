#ifndef SWAP_HPP
#define SWAP_HPP

#include "internal.hpp"

namespace Omega_h {

bool swap_part1(Mesh* mesh, AdaptOpts const& opts);

void filter_swap_improve(Mesh* mesh, LOs* cands2edges, Reals* cand_quals);

bool swap_edges(Mesh* mesh, AdaptOpts const& opts);

}  // end namespace Omega_h

#endif
