#include <Omega_h_build.hpp>
#include <Omega_h_dolfin.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 2, 2, 3);
#ifdef OMEGA_H_USE_MPI
  dolfin::Mesh mesh_dolfin(lib.world()->get_impl());
#else
  dolfin::Mesh mesh_dolfin;
#endif
  Omega_h::to_dolfin(mesh_dolfin, &mesh);
  dolfin::File file_dolfin("converted.pvd");
  file_dolfin << mesh_dolfin;
}
