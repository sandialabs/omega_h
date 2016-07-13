#include "map.hpp"
#include "mark.hpp"
#include "omega_h.hpp"
#include "surface.hpp"

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  osh::Mesh mesh;
  osh::gmsh::read(argv[1], lib, &mesh);
  auto sdim = mesh.dim() - 1;
  auto sides_are_surf = osh::mark_by_class_dim(&mesh, sdim, sdim);
  auto verts_are_surf = osh::mark_by_class_dim(&mesh, osh::VERT, sdim);
  auto surf_side2side = osh::collect_marked(sides_are_surf);
  auto surf_vert2vert = osh::collect_marked(verts_are_surf);
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
  if (mesh.dim() == 3) {
    auto edges_are_curv = osh::mark_by_class_dim(&mesh, osh::EDGE, osh::EDGE);
    auto verts_are_curv = osh::mark_by_class_dim(&mesh, osh::VERT, osh::EDGE);
    auto curv_edge2edge = osh::collect_marked(edges_are_curv);
    auto curv_vert2vert = osh::collect_marked(verts_are_curv);
    auto curv_edge_tangents =
        osh::surf::get_edge_tangents(&mesh, curv_edge2edge);
    auto edge_tangents = map_onto(curv_edge_tangents, curv_edge2edge,
                                  mesh.nedges(), 0.0, mesh.dim());
    mesh.add_tag(osh::EDGE, "tangent", mesh.dim(), OSH_DONT_TRANSFER,
                 edge_tangents);
    auto curv_vert_tangents = osh::surf::get_vert_tangents(
        &mesh, curv_edge2edge, curv_edge_tangents, curv_vert2vert);
    auto vert_tangents = map_onto(curv_vert_tangents, curv_vert2vert,
                                  mesh.nverts(), 0.0, mesh.dim());
    mesh.add_tag(osh::VERT, "tangent", mesh.dim(), OSH_DONT_TRANSFER,
                 vert_tangents);
    osh::vtk::write_vtu(argv[3], &mesh, osh::EDGE);
  }
}
