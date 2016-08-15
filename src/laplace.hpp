#ifndef LAPLACE_HPP
#define LAPLACE_HPP

#include "internal.hpp"

namespace osh {

Reals solve_laplacian(
    Mesh* mesh, Reals initial, Int width, Real tol, Real floor = EPSILON);

}  // end namespace osh

#endif
