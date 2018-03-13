#include <Omega_h_library.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_assemble.hpp>
#include <Omega_h_loop.hpp>
#include <Omega_h_shape.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_solve.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0);
  auto u_w = Omega_h::Write<Omega_h::Real>(mesh.nverts());
  auto coords = mesh.coords();
  auto f = OMEGA_H_LAMBDA(Omega_h::LO v) {
    auto x = Omega_h::get_vector<2>(coords, v);
    u_w[v] = x[0];
  };
  Omega_h::parallel_for(mesh.nverts(), f);
  auto u = Omega_h::Reals(u_w);
  auto elem_material_matrices =
    Omega_h::repeat_symm(mesh.nelems(), Omega_h::identity_matrix<2, 2>());
  auto elem_jac_invs =
    Omega_h::get_elem_jac_invs(&mesh);
  auto elem_sizes = Omega_h::measure_elements_real(&mesh);
  auto a_edge =
    Omega_h::get_edge_grad_grad(&mesh, elem_jac_invs,
        elem_material_matrices, elem_sizes);
  auto a_vert =
    Omega_h::get_vert_grad_grad(&mesh, elem_jac_invs,
        elem_material_matrices, elem_sizes);
  auto b = Omega_h::Reals(mesh.nverts(), 0.0);
  auto vert_class_dims = mesh.get_array<Omega_h::Byte>(0, "class_dim");
  auto verts_are_bdry = Omega_h::each_lt(vert_class_dims, Omega_h::Byte(2));
  mesh.add_tag(0, "a", 1, a_vert);
  mesh.add_tag(1, "a", 1, a_edge);
  mesh.add_tag(0, "u*", 1, u);
  mesh.add_tag(0, "b", 1, b);
  mesh.add_tag(0, "dbc", 1, verts_are_bdry);
  Omega_h::vtk::write_vtu("before_dbc.vtu", &mesh, 1);
  Omega_h::apply_dirichlet(&mesh, &a_edge, a_vert, &b, verts_are_bdry, u);
  mesh.set_tag(1, "a", a_edge);
  mesh.set_tag(0, "b", b);
  Omega_h::vtk::write_vtu("after_dbc.vtu", &mesh, 1);
  auto u2 = Omega_h::Reals(mesh.nverts(), 0.0);
  u2 = Omega_h::conjugate_gradient(&mesh, b, a_edge, a_vert,
      u2, 1e-6, 10);
  mesh.set_tag(0, "u", u2);
}
