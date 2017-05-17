#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_mesh.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<std::string>("output.osh");
  if (!cmdline.parse_all_or_help(world, &argc, argv)) return -1;
  auto inpath = cmdline.get<std::string>("input.osh");
  auto outpath = cmdline.get<std::string>("output.osh");
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(inpath, world, &mesh);
  Omega_h::reorder_by_hilbert(&mesh);
  Omega_h::binary::write(outpath, &mesh);
}
