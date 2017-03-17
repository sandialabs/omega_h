#include "Omega_h.hpp"
#include "Omega_h_timer.hpp"
#include "access.hpp"

#include <iostream>

using namespace Omega_h;

constexpr Int dim = 3;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh(&lib);
  auto orig_height = 4;
  auto orig_width = 1;
  auto orig_resolution = 3;
  auto max_size = 1.0 / Real(orig_resolution);
  auto segment_angle = PI / 32.0;
  if (world->rank() == 0) {
    build_box(&mesh, orig_width, orig_width, orig_height,
        orig_width * orig_resolution, orig_width * orig_resolution,
        orig_height * orig_resolution);
    classify_by_angles(&mesh, PI / 4);
    mesh.reorder();
    mesh.reset_globals();
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.add_tag(VERT, "orig_coords", mesh.dim(), OMEGA_H_LINEAR_INTERP,
      OMEGA_H_DO_OUTPUT, mesh.coords());
  mesh.add_tag<Real>(VERT, "size", 1, OMEGA_H_SIZE, OMEGA_H_DO_OUTPUT);
  vtk::Writer writer(&mesh, "bend", mesh.dim());
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
    auto f = LAMBDA(LO v) {
      auto op = get_vector<dim>(orig_coords, v);
      auto angle = full_angle * (op[2] / orig_height);
      auto radius = bend_radius - op[1];
      auto y = bend_radius - radius * cos(angle);
      auto z = radius * sin(angle);
      auto p2 = vector_3(op[0], y, z);
      auto p = get_vector<dim>(coords, v);
      auto wv = p2 - p;
      set_vector(warp_w, v, wv);
    };
    parallel_for(mesh.nverts(), f);
    mesh.add_tag(VERT, "warp", mesh.dim(), OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals(warp_w));
    do {
      std::cout << "WARP STEP\n";
      auto isos = get_curvature_isos(&mesh, segment_angle, max_size);
      isos = limit_size_field_gradation(&mesh, isos, 1.015);
      mesh.set_tag(VERT, "size", isos);
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
