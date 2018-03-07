#ifndef OMEGA_H_SOLVE_HPP
#define OMEGA_H_SOLVE_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_comm.hpp>

namespace Omega_h {

class Mesh;

Reals matrix_vector_product(Mesh* mesh, Reals a_edge, Reals a_vert, Reals x);
Real vector_dot_product(CommPtr comm, Reals u, Reals v);
Reals conjugate_gradient(Mesh* mesh, Reals b, Reals a_edge, Reals a_vert,
    Reals x_0, Real tolerance, Int max_iters);

}  // namespace Omega_h

#endif
