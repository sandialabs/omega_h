#include <Omega_h_adapt.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<double>("<target nelems>");
  cmdline.add_arg<std::string>("output.osh");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto path_in = cmdline.get<std::string>("input.osh");
  auto target_nelems = cmdline.get<double>("<target nelems>");
  auto path_out = cmdline.get<std::string>("output.osh");
  auto mesh = Omega_h::binary::read(path_in, &lib, /*strict=*/true);
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto metrics = Omega_h::get_implied_isos(&mesh);
  auto scalar =
      Omega_h::get_metric_scalar_for_nelems(&mesh, metrics, target_nelems);
  metrics = multiply_each_by(metrics, scalar);
  mesh.add_tag(Omega_h::VERT, "metric", 1, metrics);
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.verbosity = Omega_h::EXTRA_STATS;
  Omega_h::adapt(&mesh, opts);
  mesh.remove_tag(Omega_h::VERT, "metric");
  Omega_h::binary::write(path_out, &mesh);
}
