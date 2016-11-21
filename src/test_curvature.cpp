#include "Omega_h.hpp"
#include "map.hpp"
#include "mark.hpp"
#include "surface.hpp"
#include "access.hpp"
#include "space.hpp"
#include "loop.hpp"

#include <sstream>

using namespace Omega_h;

static void attach_basis_vectors(Mesh* mesh, Int ent_dim, LOs surf_ents2ents,
    Reals surf_ent_normals) {
  auto nsurf_ents = surf_ents2ents.size();
  Write<Real> surf_ent_axes[2];
  for (Int i = 0; i < 2; ++i) {
    surf_ent_axes[i] = Write<Real>(nsurf_ents * 3);
  }
  auto f = LAMBDA(LO surf_ent) {
    auto n = get_vector<3>(surf_ent_normals, surf_ent);
    auto nuv = form_ortho_basis(n);
    set_vector(surf_ent_axes[0], surf_ent, nuv[1]);
    set_vector(surf_ent_axes[1], surf_ent, nuv[2]);
  };
  parallel_for(nsurf_ents, f);
  for (Int i = 0; i < 2; ++i) {
    auto a = Reals(surf_ent_axes[i]);
    auto b = map_onto(a, surf_ents2ents, mesh->nents(ent_dim), 0.0, 3);
    std::stringstream ss;
    ss << "axis_" << i;
    auto s = ss.str();
    mesh->add_tag(ent_dim, s, 3, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT, b);
  }
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  Mesh mesh(&lib);
  CHECK(argc == 3);
  gmsh::read(argv[1], &mesh);
  auto sdim = mesh.dim() - 1;
  auto sides_are_surf = mark_by_class_dim(&mesh, sdim, sdim);
  auto verts_are_surf = mark_by_class_dim(&mesh, VERT, sdim);
  auto surf_side2side = collect_marked(sides_are_surf);
  auto surf_vert2vert = collect_marked(verts_are_surf);
  auto surf_side_normals =
      surf::get_side_normals(&mesh, surf_side2side);
  auto side_normals = map_onto(
      surf_side_normals, surf_side2side, mesh.nents(sdim), 0.0, mesh.dim());
  mesh.add_tag(sdim, "normal", mesh.dim(), OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DO_OUTPUT, side_normals);
  auto surf_vert_normals = surf::get_vert_normals(
      &mesh, surf_side2side, surf_side_normals, surf_vert2vert);
  auto vert_normals = map_onto(
      surf_vert_normals, surf_vert2vert, mesh.nverts(), 0.0, mesh.dim());
  mesh.add_tag(VERT, "normal", mesh.dim(), OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DO_OUTPUT, vert_normals);
  if (mesh.dim() == 3) {
    auto edges_are_curv =
        mark_by_class_dim(&mesh, EDGE, EDGE);
    auto verts_are_curv =
        mark_by_class_dim(&mesh, VERT, EDGE);
    auto curv_edge2edge = collect_marked(edges_are_curv);
    auto curv_vert2vert = collect_marked(verts_are_curv);
    auto curv_edge_tangents =
        surf::get_edge_tangents(&mesh, curv_edge2edge);
    auto edge_tangents = map_onto(
        curv_edge_tangents, curv_edge2edge, mesh.nedges(), 0.0, mesh.dim());
    mesh.add_tag(EDGE, "tangent", mesh.dim(), OMEGA_H_DONT_TRANSFER,
        OMEGA_H_DO_OUTPUT, edge_tangents);
    auto curv_vert_tangents = surf::get_vert_tangents(
        &mesh, curv_edge2edge, curv_edge_tangents, curv_vert2vert);
    auto vert_tangents = map_onto(
        curv_vert_tangents, curv_vert2vert, mesh.nverts(), 0.0, mesh.dim());
    mesh.add_tag(VERT, "tangent", mesh.dim(), OMEGA_H_DONT_TRANSFER,
        OMEGA_H_DO_OUTPUT, vert_tangents);
    attach_basis_vectors(&mesh, VERT, surf_vert2vert, surf_vert_normals);
    attach_basis_vectors(&mesh, TRI, surf_side2side, surf_side_normals);
    auto surf_tri_IIs = surf::get_triangle_curvatures(&mesh, surf_side2side,
        surf_side_normals, surf_vert2vert, surf_vert_normals);
    auto tri_IIs = map_onto(surf_tri_IIs, surf_side2side, mesh.ntris(),
        0.0, 3);
    mesh.add_tag(TRI, "II", 3, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT,
        tri_IIs);
    auto surf_vert_IIs = surf::get_vert_curvatures(&mesh, surf_side2side,
        surf_side_normals, surf_tri_IIs, surf_vert2vert, surf_vert_normals);
    auto vert_IIs = map_onto(surf_vert_IIs, surf_vert2vert, mesh.nverts(),
        0.0, 3);
    mesh.add_tag(VERT, "II", 3, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT,
        vert_IIs);
  }
  vtk::write_vtu(argv[2], &mesh, sdim);
}
