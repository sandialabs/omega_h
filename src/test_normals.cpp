#include "map.hpp"
#include "mark.hpp"
#include "omega_h.hpp"
#include "surface.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::Mesh mesh;
  Omega_h::gmsh::read(argv[1], lib, &mesh);
  auto sdim = mesh.dim() - 1;
  auto sides_are_surf = Omega_h::mark_by_class_dim(&mesh, sdim, sdim);
  auto verts_are_surf = Omega_h::mark_by_class_dim(&mesh, Omega_h::VERT, sdim);
  auto surf_side2side = Omega_h::collect_marked(sides_are_surf);
  auto surf_vert2vert = Omega_h::collect_marked(verts_are_surf);
  auto surf_side_normals =
      Omega_h::surf::get_side_normals(&mesh, surf_side2side);
  auto side_normals = map_onto(
      surf_side_normals, surf_side2side, mesh.nents(sdim), 0.0, mesh.dim());
  mesh.add_tag(sdim, "normal", mesh.dim(), OMEGA_H_DONT_TRANSFER, side_normals);
  auto surf_vert_normals = Omega_h::surf::get_vert_normals(
      &mesh, surf_side2side, surf_side_normals, surf_vert2vert);
  auto vert_normals = map_onto(
      surf_vert_normals, surf_vert2vert, mesh.nverts(), 0.0, mesh.dim());
  mesh.add_tag(
      Omega_h::VERT, "normal", mesh.dim(), OMEGA_H_DONT_TRANSFER, vert_normals);
  Omega_h::vtk::write_vtu(argv[2], &mesh, sdim);
  if (mesh.dim() == 3) {
    auto edges_are_curv =
        Omega_h::mark_by_class_dim(&mesh, Omega_h::EDGE, Omega_h::EDGE);
    auto verts_are_curv =
        Omega_h::mark_by_class_dim(&mesh, Omega_h::VERT, Omega_h::EDGE);
    auto curv_edge2edge = Omega_h::collect_marked(edges_are_curv);
    auto curv_vert2vert = Omega_h::collect_marked(verts_are_curv);
    auto curv_edge_tangents =
        Omega_h::surf::get_edge_tangents(&mesh, curv_edge2edge);
    auto edge_tangents = map_onto(
        curv_edge_tangents, curv_edge2edge, mesh.nedges(), 0.0, mesh.dim());
    mesh.add_tag(Omega_h::EDGE, "tangent", mesh.dim(), OMEGA_H_DONT_TRANSFER,
        edge_tangents);
    auto curv_vert_tangents = Omega_h::surf::get_vert_tangents(
        &mesh, curv_edge2edge, curv_edge_tangents, curv_vert2vert);
    auto vert_tangents = map_onto(
        curv_vert_tangents, curv_vert2vert, mesh.nverts(), 0.0, mesh.dim());
    mesh.add_tag(Omega_h::VERT, "tangent", mesh.dim(), OMEGA_H_DONT_TRANSFER,
        vert_tangents);
    Omega_h::vtk::write_vtu(argv[3], &mesh, Omega_h::EDGE);
  }
}
