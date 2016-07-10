#ifndef SWAP_HPP
#define SWAP_HPP

#include "internal.hpp"

namespace osh {

bool swap_part1(Mesh* mesh, Real qual_ceil, Int nlayers);

void filter_swap_improve(Mesh* mesh, LOs* cands2edges, Reals* cand_quals);

bool swap_edges(Mesh* mesh, Real qual_ceil, Int nlayers, bool verbose);

} //end namespace osh

#endif
