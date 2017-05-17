#ifndef OMEGA_H_LAPLACE_HPP
#define OMEGA_H_LAPLACE_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

Reals solve_laplacian(
    Mesh* mesh, Reals initial, Int width, Real tol, Real floor = EPSILON);

}  // end namespace Omega_h

#endif
