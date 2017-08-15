#include <Omega_h_adapt.hpp>
#include <Omega_h_map.hpp>
#include "Omega_h_motion.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto mesh = Mesh(&lib);
  build_box_internal(&mesh, 1.0, 1.0, 0.0, 2, 2, 0);
  mesh.add_tag(mesh.dim(), "size_error", 1, Read<Real>(mesh.nelems(), -0.1, 0.2 / 7.0));
  add_implied_metric_tag(&mesh);
  vtk::write_vtu("before.vtu", &mesh);
  auto opts = AdaptOpts(&mesh);
  opts.motion_step_size = 1.0;
  auto cands2verts = LOs({4});
  auto choices = get_motion_choices(&mesh, opts, cands2verts);
  auto new_coords = deep_copy(mesh.coords());
  map_into(choices.new_coords, cands2verts, new_coords, mesh.dim());
  mesh.set_coords(new_coords);
  vtk::write_vtu("after.vtu", &mesh);
}
