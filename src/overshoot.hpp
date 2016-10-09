#ifndef OVERSHOOT_HPP
#define OVERSHOOT_HPP

#include "Omega_h.hpp"

namespace Omega_h {

Read<I8> prevent_overshoot(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes,
    Real max_length);

}

#endif
