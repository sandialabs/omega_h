#include <Omega_h_build.hpp>
#include <Omega_h_amr.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_HYPERCUBE, 1.0, 0.0, 0.0, 1, 0, 0);
  auto xfer_opts = Omega_h::TransferOpts();
  Omega_h::amr_refine(&mesh, Omega_h::Bytes(1, Omega_h::Byte(1)), xfer_opts);
}
