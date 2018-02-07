#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_HYPERCUBE, 1.0, 1.0, 0.0, 8, 8, 0);
  Omega_h::vtk::write_vtu("test.vtu", &mesh);
}
