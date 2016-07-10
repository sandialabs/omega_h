#ifndef ADAPT_HPP
#define ADAPT_HPP

#include "internal.hpp"

namespace osh {

/* if the mesh satisfies quality and
   length requirements, print a short
   message and return true.
   otherwise, print a more detailed
   statistical report and return false */
bool adapt_check(Mesh* mesh, Real qual_floor, Real qual_ceil, Real len_floor,
                 Real len_ceil, bool verbose = true);

} //end namespace osh

#endif
