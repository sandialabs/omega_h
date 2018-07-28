#include <Omega_h_adapt.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>

#include <iostream>

using namespace Omega_h;

constexpr Int dim = 3;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  auto orig_height = 4;
  auto orig_width = 1;
  auto orig_resolution = 3;
  auto max_size = 1.0 / Real(orig_resolution);
  auto segment_angle = PI / 32.0;
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, orig_width, orig_width,
      orig_height, orig_width * orig_resolution, orig_width * orig_resolution,
      orig_height * orig_resolution);
  mesh.add_tag(VERT, "orig_coords", mesh.dim(), mesh.coords());
  mesh.add_tag<Real>(VERT, "metric", 1);
  vtk::Writer writer("bend", &mesh, mesh.dim());
  auto first_bend_radius = 5.0;
  auto final_bend_radius = orig_height / PI;
  auto nsteps = 20;
  AdaptOpts opts(&mesh);
  auto t0 = now();
  for (Int i = 0; i <= nsteps; ++i) {
    std::cout << "OUTER STEP " << i << '\n';
    auto t = Real(i) / Real(nsteps);
    auto bend_radius = first_bend_radius * (1.0 - t) + final_bend_radius * t;
    auto full_angle = orig_height / bend_radius;
    mesh.set_parting(OMEGA_H_GHOSTED);
    auto orig_coords = mesh.get_array<Real>(VERT, "orig_coords");
    auto coords = mesh.coords();
    auto warp_w = Write<Real>(mesh.nverts() * dim);
    auto f = OMEGA_H_LAMBDA(LO v) {
      auto op = get_vector<dim>(orig_coords, v);
      auto angle = full_angle * (op[2] / orig_height);
      auto radius = bend_radius - op[1];
      auto y = bend_radius - radius * std::cos(angle);
      auto z = radius * std::sin(angle);
      auto p2 = vector_3(op[0], y, z);
      auto p = get_vector<dim>(coords, v);
      auto wv = p2 - p;
      set_vector(warp_w, v, wv);
    };
    parallel_for(mesh.nverts(), f);
    mesh.add_tag(VERT, "warp", mesh.dim(), Reals(warp_w));
    do {
      std::cout << "WARP STEP\n";
      auto metric = get_curvature_metrics(&mesh, segment_angle);
      metric = apply_isotropy(mesh.nverts(), metric, OMEGA_H_ISO_LENGTH);
      metric = clamp_metrics(mesh.nverts(), metric, 0.0, max_size);
      metric = limit_metric_gradation(&mesh, metric, 1.015);
      mesh.set_tag(VERT, "metric", metric);
      adapt(&mesh, opts);
      writer.write();
    } while (warp_to_limit(&mesh, opts));
  }
  auto t1 = now();
  std::cout << "total time " << t1 - t0 << '\n';
  // paraview workaround
  writer.write();
  writer.write();
}
