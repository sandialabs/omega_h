#include <Omega_h_build.hpp>
#include <Omega_h_assemble.hpp>
#include <Omega_h_file.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1., 1., 0., 3, 4, 0);
  auto elem_material_matrices = Omega_h::repeat_symm(mesh.nelems(), Omega_h::identity_matrix<2, 2>());
  auto elem_jac_invs = Omega_h::get_elem_jac_invs(&mesh);
  auto vert_grad_grads = Omega_h::get_vert_grad_grad(&mesh, elem_jac_invs, elem_material_matrices);
  auto edge_grad_grads = Omega_h::get_edge_grad_grad(&mesh, elem_jac_invs, elem_material_matrices);
  mesh.add_tag(OMEGA_H_VERT, "grad_grad", 1, vert_grad_grads);
  mesh.add_tag(OMEGA_H_EDGE, "grad_grad", 1, edge_grad_grads);
  Omega_h::vtk::write_vtu("grad_grad.vtu", &mesh, OMEGA_H_EDGE);
}
