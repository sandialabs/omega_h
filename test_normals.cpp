#include "map.hpp"
#include "mark.hpp"
#include "omega_h.hpp"
#include "surface.hpp"

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  OSH_CHECK(argc == 3);
  osh::Mesh mesh;
  osh::gmsh::read(argv[1], lib, &mesh);
  auto sdim = mesh.dim() - 1;
  auto sides_are_surf = osh::mark_by_class_dim(&mesh, sdim, sdim);
  auto verts_are_surf = osh::mark_down(&mesh, sdim, osh::VERT, sides_are_surf);
  auto surf_side2side = collect_marked(sides_are_surf);
  auto surf_vert2vert = collect_marked(verts_are_surf);
  auto surf_side_normals = osh::surf::get_side_normals(&mesh, surf_side2side);
  auto side_normals = map_onto(surf_side_normals, surf_side2side,
                               mesh.nents(sdim), 0.0, mesh.dim());
  mesh.add_tag(sdim, "normal", mesh.dim(), OSH_DONT_TRANSFER, side_normals);
  auto surf_vert_normals = osh::surf::get_vert_normals(
      &mesh, surf_side2side, surf_side_normals, surf_vert2vert);
  auto vert_normals = map_onto(surf_vert_normals, surf_vert2vert, mesh.nverts(),
                               0.0, mesh.dim());
  mesh.add_tag(osh::VERT, "normal", mesh.dim(), OSH_DONT_TRANSFER,
               vert_normals);
  osh::vtk::write_vtu(argv[2], &mesh, sdim);
}
