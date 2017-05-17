#ifndef OMEGA_H_REFINE_HPP
#define OMEGA_H_REFINE_HPP

#include <Omega_h_adapt.hpp>

namespace Omega_h {

bool refine_by_size(Mesh* mesh, AdaptOpts const& opts);

}  // end namespace Omega_h

#endif
