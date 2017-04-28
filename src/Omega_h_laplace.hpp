#ifndef LAPLACE_HPP
#define LAPLACE_HPP


namespace Omega_h {

Reals solve_laplacian(
    Mesh* mesh, Reals initial, Int width, Real tol, Real floor = EPSILON);

}  // end namespace Omega_h

#endif
