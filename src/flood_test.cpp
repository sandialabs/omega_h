#include <Omega_h_library.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_flood.hpp>

using namespace Omega_h;

static constexpr ClassId weight_id = 20;
static constexpr ClassId droplet_id = 31;
static constexpr Int dim = 2;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 2);
  auto filename = argv[1];
  auto mesh = gmsh::read(filename, lib.world());
  mesh.set_parting(OMEGA_H_GHOSTED);
  add_implied_metric_tag(&mesh);
  auto verts_are_obj =
    mark_class_closures(&mesh, VERT, dim, {weight_id, droplet_id});
  auto verts_on_obj = collect_marked(verts_are_obj);
  auto obj_vert_warps =
    repeat_vector(verts_on_obj.size(), vector_2(0, -3./8.));
  auto vert_warps =
    map_onto(obj_vert_warps, verts_on_obj, mesh.nverts(), 0.0, dim);
  mesh.add_tag(VERT, "warp", dim, vert_warps);
  vtk::Writer writer("flood", &mesh);
  writer.write();
  auto opts = AdaptOpts(&mesh);
  while (warp_to_limit(&mesh, opts)) {
    adapt(&mesh, opts);
    if (mesh.min_quality() < opts.min_quality_desired) {
      writer.write();
      return 0;
    }
  }
}
