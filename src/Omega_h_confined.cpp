#include "Omega_h_confined.hpp"

#include "access.hpp"
#include "internal.hpp"

/* Code to find classification-enforced constraints on
 * element size and shape.
 * For example, two geometric boundaries may be so close
 * as to force a mesh edge to be much smaller than desired,
 * similarly a sharp angle may be so acute that elements
 * contained by it are forced to be low-quality
 */

namespace Omega_h {

Read<I8> find_bridge_edges(Mesh* mesh) {
  auto verts2class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2class_dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto edges_are_bridges_w = Write<I8>(mesh->nedges());
  auto edges2verts = mesh->ask_verts_of(EDGE);
  auto f = LAMBDA(LO edge) {
    auto eev2v = gather_verts<2>(edges2verts, edge);
    auto eev2dim = gather_scalars(edges2class_dim, eev2v);
    auto edim = edges2class_dim[edge];
    edges_are_bridges_w[edge] = ((edim != verts2class_dim[eev2v[0]]) &&
                                 (edim != verts2class_dim[eev2v[1]]));
  };
  parallel_for(mesh->nedges(), f);
  return edges_are_bridges_w;
}

static DEVICE bool is_angle_triangle(
    Few<I8, 3> vert_dims, Few<I8, 3> edge_dims, I8 tri_dim) {
  for (Int i = 0; i < 3; ++i) {
    if (vert_dims[i] > (tri_dim - 2)) continue;
    if (edge_dims[(i + 0) % 3] > (tri_dim - 1)) continue;
    if (edge_dims[(i + 2) % 3] > (tri_dim - 1)) continue;
    return true;
  }
  return false;
}

Read<I8> find_angle_triangles(Mesh* mesh) {
  auto verts2class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2class_dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto tris2class_dim = mesh->get_array<I8>(TRI, "class_dim");
  auto tris2edges = mesh->ask_down(TRI, EDGE).ab2b;
  auto tris2verts = mesh->ask_down(TRI, VERT).ab2b;
  auto tris_are_angle = Write<I8>(mesh->ntris());
  auto f = LAMBDA(LO tri) {
    auto ttv2v = gather_down<3>(tris2verts, tri);
    auto tte2e = gather_down<3>(tris2edges, tri);
    auto ttv2dim = gather_scalars(verts2class_dim, ttv2v);
    auto tte2dim = gather_scalars(edges2class_dim, tte2e);
    auto t_dim = tris2class_dim[tri];
    tris_are_angle[tri] = is_angle_triangle(ttv2dim, tte2dim, t_dim);
  };
  parallel_for(mesh->ntris(), f);
  return tris_are_angle;
}

Read<I8> find_angle_elems(Mesh* mesh) {
  auto tris_are_angle = find_angle_triangles(mesh);
  if (mesh->dim() == 2) return tris_are_angle;
  if (mesh->dim() == 3) return mark_up(mesh, TRI, TET, tris_are_angle);
  NORETURN(Read<I8>());
}

}  // end namespace Omega_h
