#include "omega_h.hpp"

#include "access.hpp"
#include "internal.hpp"
#include "loop.hpp"
#include "derive.hpp"
#include "project.hpp"

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  osh::Mesh mesh;
  osh::build_box(&mesh, lib, 1, 1, 0, 8, 8, 0);
  osh::classify_by_angles(&mesh, osh::PI / 4);
  auto u_w = osh::Write<osh::Real>(mesh.nverts());
  auto coords = mesh.coords();
  auto f = LAMBDA(osh::LO v) {
    auto x = osh::get_vector<2>(coords, v);
    u_w[v] = osh::square(x[0]) + osh::square(x[1]);
  };
  osh::parallel_for(mesh.nverts(), f);
  auto u = osh::Reals(u_w);
  mesh.add_tag(osh::VERT, "u", 1, OSH_DONT_TRANSFER, u);
  auto e_grad = osh::derive_element_gradients(&mesh, u);
  mesh.add_tag(mesh.dim(), "e_grad", mesh.dim(), OSH_DONT_TRANSFER, e_grad);
  auto grad = osh::project(&mesh, e_grad);
  mesh.add_tag(osh::VERT, "grad", mesh.dim(), OSH_DONT_TRANSFER, grad);
  auto e_hess = osh::derive_element_hessians(&mesh, grad);
  mesh.add_tag(mesh.dim(), "e_hess", osh::symm_dofs(mesh.dim()), OSH_DONT_TRANSFER, e_hess);
  auto hess = osh::project(&mesh, e_hess);
  mesh.add_tag(osh::VERT, "hess", osh::symm_dofs(mesh.dim()), OSH_DONT_TRANSFER, hess);
  osh::vtk::write_vtu("hessian.vtu", &mesh, mesh.dim());
}
