#ifndef REFINE_HPP
#define REFINE_HPP

#include "internal.hpp"

namespace Omega_h {

/* This is sqrt(2), but sqrt() is not constexpr */
constexpr Real max_length_desired = 1.4142135623730951;

bool refine_by_size(Mesh* mesh, AdaptOpts const& opts);

}  // end namespace Omega_h

#endif
