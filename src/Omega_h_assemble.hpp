#ifndef OMEGA_H_ASSEMBLE_HPP
#define OMEGA_H_ASSEMBLE_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

Reals get_edge_grad_grad(Mesh* mesh, Reals elem_basis_grads,
    Reals elem_material_matrices, Reals elem_sizes);
Reals get_vert_grad_grad(Mesh* mesh, Reals elem_basis_grads,
    Reals elem_material_matrices, Reals elem_sizes);
Reals get_elem_jac_invs(Mesh* mesh);
void apply_dirichlet(Mesh* mesh, Reals* p_a_edge, Reals a_vert, Reals* p_b,
    Bytes are_dbc, Reals dbc_values);

}  // namespace Omega_h

#endif
