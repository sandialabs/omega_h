#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.osh");
  cmdline.add_arg<std::string>("mesh_out.osh");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto path_in = cmdline.get<std::string>("mesh_in.osh");
  auto path_out = cmdline.get<std::string>("mesh_out.osh");
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(path_in, lib.world(), &mesh);
  mesh.remove_tag(Omega_h::VERT, "metric");
  Omega_h::add_implied_metric_tag(&mesh);
  mesh.ask_lengths();
  mesh.ask_qualities();
  Omega_h::binary::write(path_out, &mesh);
}
