#include <iostream>

#include "Omega_h_adapt.hpp"
#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_metric.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
#ifdef OMEGA_H_USE_LIBMESHB
  auto metric_doc = "metric.{txt,sol}";
#else
  auto metric_doc = "metric.txt";
#endif
  auto& mesh_in_flag = cmdline.add_flag("--mesh-in", "input mesh file path");
  mesh_in_flag.add_arg<std::string>("mesh.msh");
  auto& mesh_out_flag = cmdline.add_flag("--mesh-out", "output mesh file path");
  mesh_out_flag.add_arg<std::string>("mesh.msh");
  auto& metric_in_flag =
      cmdline.add_flag("--metric-in", "input metric file path");
  metric_in_flag.add_arg<std::string>(metric_doc);
  auto& metric_out_flag =
      cmdline.add_flag("--metric-out", "output metric file path");
  metric_out_flag.add_arg<std::string>(metric_doc);
  if (!cmdline.parse_final(comm, &argc, argv)) {
    return -1;
  }
  if (!cmdline.parsed("--mesh-in")) {
    std::cout << "No input mesh file specified\n";
    cmdline.show_help(argv);
    return -1;
  }
  if (!cmdline.parsed("--mesh-out")) {
    std::cout << "No output mesh file specified\n";
    cmdline.show_help(argv);
    return -1;
  }
  if (!cmdline.parsed("--metric-in")) {
    std::cout << "No input metric file specified\n";
    cmdline.show_help(argv);
    return -1;
  }
  Omega_h::filesystem::path mesh_in =
      cmdline.get<std::string>("--mesh-in", "mesh.msh");
  Omega_h::filesystem::path mesh_out =
      cmdline.get<std::string>("--mesh-out", "mesh.msh");
  Omega_h::filesystem::path metric_in =
      cmdline.get<std::string>("--metric-in", metric_doc);
  std::cout << "Loading mesh from " << mesh_in << "\n";
  auto mesh = Omega_h::gmsh::read(mesh_in, comm);
  Omega_h::Reals target_metric;
  std::cout << "Loading target metric from " << metric_in << "\n";
  auto const dim = mesh.dim();
  auto const metric_in_ext = metric_in.extension().string();
  if (metric_in_ext == ".txt") {
    target_metric = Omega_h::read_reals_txt(
        metric_in, mesh.nverts(), Omega_h::symm_ncomps(dim));
    target_metric = Omega_h::symms_inria2osh(dim, target_metric);
    mesh.add_tag(0, "target_metric", Omega_h::symm_ncomps(dim), target_metric);
  } else
#ifdef OMEGA_H_USE_LIBMESHB
      if (metric_in_ext == ".sol" || metric_in_ext == ".solb") {
    Omega_h::meshb::read_sol(&mesh, metric_in.c_str(), "target_metric");
    target_metric = mesh.get_array<Omega_h::Real>(0, "target_metric");
  } else
#endif
  {
    Omega_h_fail("unknown extension for \"%s\"\n", metric_in.c_str());
  }
  auto opts = Omega_h::AdaptOpts(&mesh);
  Omega_h::grade_fix_adapt(&mesh, opts, target_metric, /*verbose=*/true);
  std::cout << "Storing mesh in " << mesh_out << '\n';
  Omega_h::gmsh::write(mesh_out, &mesh);
  if (cmdline.parsed("--metric-out")) {
    Omega_h::filesystem::path metric_out =
        cmdline.get<std::string>("--metric-out", metric_doc);
    auto const metric_out_ext = metric_out.extension().string();
    std::cout << "Storing metric in " << metric_out << '\n';
    if (metric_out_ext == ".txt") {
      auto metric = mesh.get_array<Omega_h::Real>(0, "metric");
      metric = Omega_h::symms_osh2inria(dim, metric);
      Omega_h::write_reals_txt(metric_out, metric, Omega_h::symm_ncomps(dim));
    } else
#ifdef OMEGA_H_USE_LIBMESHB
        if (metric_out_ext == ".sol" || metric_out_ext == ".solb") {
      Omega_h::meshb::write_sol(&mesh, metric_out.c_str(), "metric");
    } else
#endif
    {
      Omega_h_fail("unknown extension for \"%s\"\n", metric_out.c_str());
    }
  }
}
