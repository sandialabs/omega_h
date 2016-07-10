#ifndef REFINE_HPP
#define REFINE_HPP

#include "internal.hpp"

namespace osh {

bool refine(Mesh* mesh, Real min_qual, bool verbose);
bool refine_by_size(Mesh* mesh, Real max_len, Real min_qual, bool verbose);

} //end namespace osh

#endif
