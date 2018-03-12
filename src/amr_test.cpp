#include <Omega_h_build.hpp>
#include <Omega_h_amr.hpp>
#include <Omega_h_print.hpp>
#include <Omega_h_array_ops.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_HYPERCUBE, 1.0, 0.0, 0.0, 1, 0, 0);
  auto xfer_opts = Omega_h::TransferOpts();
  Omega_h::amr_refine(&mesh, Omega_h::Bytes(1, Omega_h::Byte(1)), xfer_opts);
  std::cout << "after refine, " << mesh.nelems() << " elements\n";
  OMEGA_H_CHECK(mesh.nelems() == 3);
  std::cout << "elems2verts: " << mesh.ask_elem_verts() << '\n';
  OMEGA_H_CHECK(mesh.ask_elem_verts() == Omega_h::LOs({0, 2, 0, 1, 1, 2}));
  std::cout << "coords: " << mesh.coords() << '\n';
  OMEGA_H_CHECK(mesh.coords() == Omega_h::Reals({0.0, 0.5, 1.0}));
}
