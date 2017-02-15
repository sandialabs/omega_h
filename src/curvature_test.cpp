#include "Omega_h.hpp"
#include "Omega_h_compare.hpp"
#include "access.hpp"
#include "eigen.hpp"
#include "host_few.hpp"
#include "loop.hpp"
#include "map.hpp"
#include "mark.hpp"
#include "space.hpp"
#include "surface.hpp"

#include <sstream>

using namespace Omega_h;

static void attach_basis_vectors(
    Mesh* mesh, Int ent_dim, LOs surf_ents2ents, Reals surf_ent_normals) {
  auto nsurf_ents = surf_ents2ents.size();
  HostFew<Write<Real>, 2> surf_ent_axes;
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
  CHECK(argc == 3);
  std::string path = argv[1];
  std::string name = argv[2];
  Mesh mesh(&lib);
  auto world = lib.world();
  if (world->rank() == 0) {
    gmsh::read(path + "/" + name + ".msh", &mesh);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto sdim = mesh.dim() - 1;
  auto sides_are_surf = mark_by_class_dim(&mesh, sdim, sdim);
  auto verts_are_surf = mark_by_class_dim(&mesh, VERT, sdim);
  auto surf_side2side = collect_marked(sides_are_surf);
  auto surf_vert2vert = collect_marked(verts_are_surf);
  auto vert_curvatures_w = Write<Real>(mesh.nverts(), 0.0);
  LOs curv_edge2edge;
  LOs curv_vert2vert;
  if (mesh.dim() == 3) {
    auto surf_side_normals = get_side_normals(&mesh, surf_side2side);
    auto side_normals = map_onto(
        surf_side_normals, surf_side2side, mesh.nents(sdim), 0.0, mesh.dim());
    mesh.add_tag(sdim, "normal", mesh.dim(), OMEGA_H_DONT_TRANSFER,
        OMEGA_H_DO_OUTPUT, side_normals);
    auto surf_vert_normals = get_vert_normals(
        &mesh, surf_side2side, surf_side_normals, surf_vert2vert);
    auto vert_normals = map_onto(
        surf_vert_normals, surf_vert2vert, mesh.nverts(), 0.0, mesh.dim());
    mesh.add_tag(VERT, "normal", mesh.dim(), OMEGA_H_DONT_TRANSFER,
        OMEGA_H_DO_OUTPUT, vert_normals);
    auto edges_are_curv = mark_by_class_dim(&mesh, EDGE, EDGE);
    auto verts_are_curv = mark_by_class_dim(&mesh, VERT, EDGE);
    curv_edge2edge = collect_marked(edges_are_curv);
    curv_vert2vert = collect_marked(verts_are_curv);
    attach_basis_vectors(&mesh, VERT, surf_vert2vert, surf_vert_normals);
    attach_basis_vectors(&mesh, TRI, surf_side2side, surf_side_normals);
    auto surf_tri_IIs = get_triangle_IIs(&mesh, surf_side2side,
        surf_side_normals, surf_vert2vert, surf_vert_normals);
    auto tri_IIs = map_onto(surf_tri_IIs, surf_side2side, mesh.ntris(), 0.0, 3);
    mesh.add_tag(
        TRI, "II", 3, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT, tri_IIs);
    auto surf_vert_IIs = get_vert_IIs(&mesh, surf_side2side, surf_side_normals,
        surf_tri_IIs, surf_vert2vert, surf_vert_normals);
    auto vert_IIs =
        map_onto(surf_vert_IIs, surf_vert2vert, mesh.nverts(), 0.0, 3);
    mesh.add_tag(
        VERT, "II", 3, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT, vert_IIs);
    auto surf_vert_curvatures = get_max_eigenvalues(2, surf_vert_IIs);
    map_into(surf_vert_curvatures, surf_vert2vert, vert_curvatures_w, 1);
  } else {
    curv_edge2edge = surf_side2side;
    curv_vert2vert = surf_vert2vert;
  }
  auto curv_edge_tangents = get_edge_tangents(&mesh, curv_edge2edge);
  auto edge_tangents = map_onto(
      curv_edge_tangents, curv_edge2edge, mesh.nedges(), 0.0, mesh.dim());
  mesh.add_tag(EDGE, "tangent", mesh.dim(), OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DO_OUTPUT, edge_tangents);
  auto curv_vert_tangents = get_vert_tangents(
      &mesh, curv_edge2edge, curv_edge_tangents, curv_vert2vert);
  auto vert_tangents = map_onto(
      curv_vert_tangents, curv_vert2vert, mesh.nverts(), 0.0, mesh.dim());
  mesh.add_tag(VERT, "tangent", mesh.dim(), OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DO_OUTPUT, vert_tangents);
  auto curv_edge_curvatures = get_edge_curvatures(&mesh, curv_edge2edge,
      curv_edge_tangents, curv_vert2vert, curv_vert_tangents);
  auto edge_curvatures =
      map_onto(curv_edge_curvatures, curv_edge2edge, mesh.nedges(), 0.0, 1);
  mesh.add_tag(EDGE, "curvature", 1, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT,
      edge_curvatures);
  auto curv_vert_curvatures = get_curv_vert_curvatures(
      &mesh, curv_edge2edge, curv_edge_curvatures, curv_vert2vert);
  map_into(curv_vert_curvatures, curv_vert2vert, vert_curvatures_w, 1);
  auto vert_curvatures = get_corner_vert_curvatures(&mesh, vert_curvatures_w);
  mesh.add_tag(VERT, "curvature", 1, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT,
      vert_curvatures);
  bool ok = check_regression(std::string("gold_curv_") + name, &mesh);
  if (!ok) return 2;
  return 0;
}
