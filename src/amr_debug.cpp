#include <Omega_h_amr.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>

using Omega_h::Bytes;

static void test_2D_arrays(Omega_h::Library* lib) {
  auto w = lib->world();
  auto f = OMEGA_H_HYPERCUBE;
  auto m = Omega_h::build_box(w, f, 2.0, 1.0, 0.0, 2, 1, 0);
  Omega_h::vtk::FullWriter writer("out", &m);
  writer.write();
  auto xfer_opts = Omega_h::TransferOpts();
  Omega_h::amr_refine(&m, Omega_h::Bytes({1, 0}), xfer_opts);
  writer.write();
  Omega_h::amr_refine(&m, Omega_h::Bytes({0,0,0,0,0,1}), xfer_opts);
  writer.write();
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  test_2D_arrays(&lib);
}
