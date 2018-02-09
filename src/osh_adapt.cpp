#include <iostream>

#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_adapt.hpp"
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
  auto& metric_in_flag = cmdline.add_flag("--metric-in", "input metric file path");
  metric_in_flag.add_arg<std::string>(metric_doc);
  auto& metric_out_flag = cmdline.add_flag("--metric-out", "output metric file path");
  metric_out_flag.add_arg<std::string>(metric_doc);
  if (!cmdline.parse_final(comm, &argc, argv)) {
    return -1;
  }
  if (!cmdline.parsed("--mesh-in")) {
    Omega_h_fail("No input mesh file specified");
  }
  if (!cmdline.parsed("--mesh-out")) {
    Omega_h_fail("No output mesh file specified");
  }
  if (!cmdline.parsed("--metric-in")) {
    Omega_h_fail("No input metric file specified");
  }
  auto mesh_in = cmdline.get<std::string>("--mesh-in", "mesh.msh");
  auto mesh_out = cmdline.get<std::string>("--mesh-out", "mesh.msh");
  auto metric_in = cmdline.get<std::string>("--metric-in", metric_doc);
  std::cout << "Loading mesh from " << mesh_in << "\n";
  auto mesh = Omega_h::gmsh::read(mesh_in, comm);
  Omega_h::Reals target_metric;
  std::cout << "Loading metric from " << metric_in << "\n";
  if (Omega_h::ends_with(metric_in, ".txt")) {
  //target_metric = Omega_h::read_reals_txt(metric_in, mesh.nverts(), Omega_h::symm_ndofs(mesh.dim()));
    target_metric = Omega_h::symms_inria2osh(mesh.dim(), target_metric);
    mesh.add_tag(0, "target_metric", Omega_h::symm_ncomps(mesh.dim()), target_metric);
  } else
#ifdef OMEGA_H_USE_LIBMESHB
  if (Omega_h::ends_with(metric_in, ".sol") ||
      Omega_h::ends_with(metric_in, ".solb")) {
    Omega_h::meshb::read_sol(&mesh, metric_in, "target_metric");
  } else
#endif
  {
    Omega_h_fail("unknown extension for \"%s\"\n", metric_in.c_str());
  }
  target_metric = Omega_h::limit_metric_gradation(&mesh, target_metric, 1.0);
  auto opts = Omega_h::AdaptOpts(&mesh);
  auto min_qual = mesh.min_quality();
  if (min_qual < opts.min_quality_allowed) {
    std::cout << "Initial mesh has minimum quality " << min_qual
      << " < minimum acceptable quality " << opts.min_quality_allowed << '\n';
    std::cout << "Omega_h will now attempt to repair the initial mesh quality.\n";
    std::cout << "This could take some time...\n";
    Omega_h::fix(&mesh, opts, OMEGA_H_ANISOTROPIC, /*verbose=*/true);
    std::cout << "\nOmega_h is done repairing mesh quality!\n\n";
  } else {
    Omega_h::add_implied_metric_tag(&mesh);
  }
  std::cout << "Adapting...\n";
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
  }
  std::cout << "\nDone adapting!\n\n";
  std::cout << "Storing mesh in " << mesh_out << '\n';
  Omega_h::gmsh::write(mesh_out, &mesh);
  if (cmdline.parsed("--metric-out")) {
    auto metric_out = cmdline.get<std::string>("--metric-out", metric_doc);
    std::cout << "Storing metric in " << metric_out << '\n';
    if (Omega_h::ends_with(metric_out, ".txt")) {
      auto metric = mesh.get_array<Omega_h::Real>(0, "metric");
      metric = Omega_h::symms_osh2inria(mesh.dim(), metric);
      //Omega_h::write_reals_txt(metric_out, metric);
    } else
#ifdef OMEGA_H_USE_LIBMESHB
    if (Omega_h::ends_with(metric_out, ".sol") ||
        Omega_h::ends_with(metric_out, ".solb")) {
      Omega_h::meshb::write_sol(&mesh, metric_out, "metric");
    } else
#endif
    {
      Omega_h_fail("unknown extension for \"%s\"\n", metric_out.c_str());
    }
  }
}
