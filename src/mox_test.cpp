#include <Omega_h_adapt.hpp>
#include <Omega_h_adj.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_recover.hpp>

#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  if (argc == 2) {
    Omega_h::binary::read(argv[1], world, &mesh);
  } else {
    mesh =
        Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 12, 12, 12);
  }
  mesh.balance();
  Omega_h::vtk::Writer writer("adapting", &mesh);
  auto tau_start = 10.0;
  auto tau_end = 1.0;
  auto niter = 30;
  for (int iter = 0; iter < niter; ++iter) {
    mesh.set_parting(OMEGA_H_GHOSTED);
    auto u_w = Omega_h::Write<Omega_h::Real>(mesh.nverts());
    auto eps = 0.01;
    auto coords = mesh.coords();
    auto u_f = OMEGA_H_LAMBDA(Omega_h::LO v) {
      auto x = coords[v * 3 + 0];
      auto y = coords[v * 3 + 1];
      auto z = coords[v * 3 + 2];
      u_w[v] = std::exp(-x / eps) + std::exp(-y / eps) + std::exp(-z / eps);
    };
    Omega_h::parallel_for(mesh.nverts(), u_f);
    auto u = Omega_h::Reals(u_w);
    mesh.remove_tag(0, "u");
    mesh.add_tag(0, "u", 1, u);
    auto gradients = Omega_h::derive_element_gradients(&mesh, u);
    auto pseudo_time = Omega_h::Real(iter) / Omega_h::Real(niter - 1);
    auto tau = tau_start * (1.0 - pseudo_time) + tau_end * pseudo_time;
    auto max_size = 1.0;
    std::cout << "using tau = " << tau << '\n';
    auto metrics =
        Omega_h::get_aniso_zz_metric(&mesh, gradients, tau, max_size);
    mesh.remove_tag(0, "original_metric");
    mesh.add_tag(0, "original_metric", 6, metrics);
    auto metrics2 = Omega_h::limit_metric_gradation(&mesh, metrics, 0.6);
    mesh.remove_tag(0, "target_metric");
    mesh.add_tag(0, "target_metric", 6, metrics2);
    Omega_h::vtk::write_vtu("metric.vtu", &mesh);
    if (!mesh.has_tag(0, "metric")) Omega_h::add_implied_metric_tag(&mesh);
    auto opts = Omega_h::AdaptOpts(&mesh);
    opts.max_length_allowed = opts.max_length_desired * 2;
    opts.verbosity = Omega_h::EXTRA_STATS;
    opts.xfer_opts.type_map["original_metric"] = OMEGA_H_METRIC;
    opts.xfer_opts.type_map["u"] = OMEGA_H_LINEAR_INTERP;
    while (Omega_h::approach_metric(&mesh, opts)) {
      Omega_h::adapt(&mesh, opts);
      writer.write();
    }
    std::cout << "END OF ITER " << iter << '\n';
    Omega_h::binary::write("restart.osh", &mesh);
    std::cout << "WROTE restart.osh\n";
  }
}
