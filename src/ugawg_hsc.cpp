#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_adapt.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, "hsc01.meshb");
  Omega_h::meshb::read_sol(&mesh, "hsc01-mach.solb", "mach");
  Omega_h::meshb::read_sol(&mesh, "hsc01-metric.solb", "target_metric");
  Omega_h::add_implied_metric_tag(&mesh);
  mesh.ask_qualities();
  Omega_h::vtk::write_vtu("hsc_preview.vtu", &mesh);
  return 0;
}

