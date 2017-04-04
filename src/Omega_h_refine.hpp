#ifndef REFINE_HPP
#define REFINE_HPP

#include "Omega_h_internal.hpp"

namespace Omega_h {

bool refine_by_size(Mesh* mesh, AdaptOpts const& opts);

}  // end namespace Omega_h

#endif
