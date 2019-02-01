#ifndef OMEGA_H_RANDOM_HPP
#define OMEGA_H_RANDOM_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

Reals unit_uniform_random_reals_from_globals(
    GOs const globals, I64 const seed, I64 const counter);

}

#endif
