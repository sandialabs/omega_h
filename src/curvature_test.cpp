#include "Omega_h.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_eigen.hpp"
#include "Omega_h_host_few.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_surface.hpp"

#include <sstream>

using namespace Omega_h;

static void attach_basis_vectors(
    Mesh* mesh, Int ent_dim, LOs surf_ents2ents, Reals surf_ent_normals) {
  auto nsurf_ents = surf_ents2ents.size();
  HostFew<Write<Real>, 2> surf_ent_axes;
  for (Int i = 0; i < 2; ++i) {
    surf_ent_axes[i] = Write<Real>(nsurf_ents * 3);
  }
  auto f = OMEGA_H_LAMBDA(LO surf_ent) {
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
    mesh->add_tag(ent_dim, s, 3, b);
  }
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 3);
  std::string path = argv[1];
  std::string name = argv[2];
  auto world = lib.world();
  auto mesh = gmsh::read(path + "/" + name + ".msh", world);
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto surface_info = get_surface_info(&mesh);
  auto vert_normals = map_onto(surface_info.surf_vert_normals,
      surface_info.surf_vert2vert, mesh.nverts(), 0.0, mesh.dim());
  mesh.add_tag(VERT, "normal", mesh.dim(), vert_normals);
  if (mesh.dim() == 3) {
    auto vert_IIs =
        map_onto(surface_info.surf_vert_IIs,
           surface_info.surf_vert2vert, mesh.nverts(), 0.0, 3);
    mesh.add_tag(VERT, "II", 3, vert_IIs);
    attach_basis_vectors(&mesh, 2, surface_info.surf_vert2vert,
        surface_info.surf_vert_normals);
  }
  auto vert_tangents = map_onto(surface_info.curv_vert_tangents,
      surface_info.curv_vert2vert, mesh.nverts(), 0.0, mesh.dim());
  mesh.add_tag(VERT, "tangent", mesh.dim(), vert_tangents);
  auto vert_curvatures = get_vert_curvatures(&mesh, surface_info);
  mesh.add_tag(VERT, "curvature", 1, vert_curvatures);
  bool ok = check_regression(std::string("gold_curv_") + name, &mesh);
  if (!ok) return 2;
  return 0;
}
