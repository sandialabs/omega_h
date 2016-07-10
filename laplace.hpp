#ifndef LAPLACE_HPP
#define LAPLACE_HPP

namespace osh {

Reals solve_laplacian(Mesh* mesh, Reals initial, Int width, Real tol);

} //end namespace osh

#endif
