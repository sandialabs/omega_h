#ifndef SIMPLIFY_HPP
#define SIMPLIFY_HPP

#include "internal.hpp"

namespace osh {

namespace simplify {
LOs tris_from_quads(LOs qv2v);
LOs tets_from_hexes(LOs hv2v);
}

}  // end namespace osh

#endif
