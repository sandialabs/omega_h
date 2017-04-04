#ifndef SIMPLIFY_HPP
#define SIMPLIFY_HPP

#include "Omega_h_internal.hpp"

namespace Omega_h {

namespace simplify {
LOs tris_from_quads(LOs qv2v);
LOs tets_from_hexes(LOs hv2v);
}  // namespace simplify

}  // end namespace Omega_h

#endif
