#ifndef OMEGA_H_COMPARE_HPP
#define OMEGA_H_COMPARE_HPP

#include "Omega_h.hpp"

namespace Omega_h {

Omega_h_Comparison compare_meshes(
    Mesh* a, Mesh* b, Real tol, Real floor, bool verbose, bool full = true);

}

#endif
