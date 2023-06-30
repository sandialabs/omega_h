#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.exo");
  cmdline.add_arg<std::string>("output.osh");
  cmdline.add_flag("-v", "verbose");
  auto& cwflag = cmdline.add_flag(
      "--classify-with", "comma separated \"node_sets,side_sets\"");
  cwflag.add_arg<std::string>("set-types");
  if (!cmdline.parse(comm, &argc, argv) ||
      !Omega_h::CmdLine::check_empty(comm, argc, argv)) {
    cmdline.show_help(comm, argv);
    return -1;
  }
  auto inpath = cmdline.get<std::string>("input.exo");
  auto outpath = cmdline.get<std::string>("output.osh");
  auto verbose = cmdline.parsed("-v");
  int classify_with;
  if (cmdline.parsed("--classify-with")) {
    auto set_types = cmdline.get<std::string>("--classify-with", "set-types");
    classify_with = 0;
    if (set_types.find("node_sets") != std::string::npos) {
      classify_with |= Omega_h::exodus::NODE_SETS;
    }
    if (set_types.find("side_sets") != std::string::npos) {
      classify_with |= Omega_h::exodus::SIDE_SETS;
    }
  } else {
    classify_with = Omega_h::exodus::NODE_SETS | Omega_h::exodus::SIDE_SETS;
  }
  Omega_h::Mesh mesh(&lib);
  mesh.set_comm(comm);
  if (comm->size() == 1) {
    auto exodus_file = Omega_h::exodus::open(inpath, verbose);
    Omega_h::exodus::read_mesh(exodus_file, &mesh, verbose, classify_with);
    auto ntime_steps = Omega_h::exodus::get_num_time_steps(exodus_file);
    if (ntime_steps > 0) {
      Omega_h::exodus::read_nodal_fields(
          exodus_file, &mesh, ntime_steps - 1, "", "", verbose);
    }
    Omega_h::exodus::close(exodus_file);
  } else {
    auto time_step = -1;
    mesh = Omega_h::exodus::read_sliced(
        inpath, comm, verbose, classify_with, time_step);
  }
  Omega_h::binary::write(outpath, &mesh);
  return 0;
}
