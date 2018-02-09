#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
#ifdef OMEGA_H_USE_LIBMESHB
  auto metricfile_placeholder = "metric.{txt,sol}";
#else
  auto metricfile_placeholder = "metric.txt";
#endif
  auto& meshfile_flag = cmdline.add_flag("--mesh-file");
  meshfile_flag.add_arg<std::string>("mesh.msh");
  auto& metricfile_flag = cmdline.add_flag("--metric-file");
  metricfile_flag.add_arg<std::string>(metricfile_placeholder);
  if (!cmdline.parse_final(comm, &argc, argv)) {
    return -1;
  }
  if (!cmdline.parsed("--mesh-file")) {
    Omega_h_fail("No mesh file specified");
  }
  if (!cmdline.parsed("--metric-file")) {
    Omega_h_fail("No metric file specified");
  }
  auto meshfile = cmdline.get<std::string>("--mesh-file", "mesh.msh");
  auto metricfile = cmdline.get<std::string>("--metric-file", metricfile_placeholder);
}
