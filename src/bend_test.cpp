#include "Omega_h.hpp"
#include "access.hpp"
#include "loop.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh(&lib);
  auto orig_height = 4;
  auto orig_width = 1;
  auto orig_resolution = 3;
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
  for (Int i = 0; i <= nsteps; ++i) {
    auto t = Real(i) / Real(nsteps);
    auto bend_radius = first_bend_radius * (1.0 - t) + final_bend_radius * t;
    auto full_angle = orig_height / bend_radius;
    mesh.set_parting(OMEGA_H_GHOSTED);
    auto orig_coords = mesh.get_array<Real>(VERT, "orig_coords");
    auto new_coords_w = Write<Real>(mesh.nverts() * 3);
    auto f = LAMBDA(LO v) {
      auto p = get_vector<3>(orig_coords, v);
      auto angle = full_angle * (p[2] / orig_height);
      auto radius = bend_radius - p[1];
      auto y = bend_radius - radius * cos(angle);
      auto z = radius * sin(angle);
      auto p2 = vector_3(p[0], y, z);
      set_vector(new_coords_w, v, p2);
    };
    parallel_for(mesh.nverts(), f);
    mesh.set_coords(Reals(new_coords_w));
    auto max_size = 1.0 / Real(orig_resolution);
    auto segment_angle = PI / 32.0;
    auto isos = get_curvature_isos(&mesh, segment_angle, max_size);
    isos = limit_size_field_gradation(&mesh, isos, 1.05);
    mesh.set_tag(VERT, "size", isos);
    writer.write();
  }
}
