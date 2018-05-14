#include <Omega_h_amr.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>

using Omega_h::Bytes;

static void test_3D_refine(Omega_h::Library* lib) {
  auto w = lib->world();
  auto f = OMEGA_H_HYPERCUBE;
  auto m = Omega_h::build_box(w, f, 1.0, 1.0, 1.0, 2, 1, 1);
  Omega_h::vtk::FullWriter writer("out", &m);
  writer.write();
  auto xfer_opts = Omega_h::TransferOpts();

  std::cout << "refine 1: " << std::endl;
  auto mark1 = Omega_h::Bytes({1, 0});
  auto one_level_mark1 = Omega_h::enforce_one_level(&m, 2, mark1);
  Omega_h::amr_refine(&m, one_level_mark1, xfer_opts);
  writer.write();

  std::cout << "refine 2: " << std::endl;
  auto mark2 = Omega_h::Bytes({0, 0, 1, 0, 0, 0, 0, 0, 0, 0});
  auto one_level_mark2 = Omega_h::enforce_one_level(&m, 2, mark2);
  Omega_h::amr_refine(&m, one_level_mark2, xfer_opts);
  writer.write();
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  test_3D_refine(&lib);
}
