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
    build_box(&mesh,
        orig_width,
        orig_width,
        orig_height,
        orig_width * orig_resolution,
        orig_width * orig_resolution,
        orig_height * orig_resolution);
    classify_by_angles(&mesh, PI / 4);
    mesh.reorder();
    mesh.reset_globals();
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto orig_coords = mesh.coords();
  auto bend_radius = 10.0;
  auto full_angle = orig_height / bend_radius;
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
  vtk::write_vtu("debug.vtu", &mesh, mesh.dim());
}
