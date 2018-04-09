#include "Omega_h_simplify.hpp"

#include "Omega_h_mesh.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_adj.hpp"
#include "Omega_h_hypercube.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_build.hpp"

namespace Omega_h {

void quads2tris(Mesh* mesh) {
  auto nold_verts = mesh->nverts();
  auto nquads = mesh->nfaces();
  auto old_coords = mesh->coords();
  auto quad_center_coords = average_field(mesh, 2, 2, old_coords);
  auto nverts = nold_verts + nquads;
  auto coords = Write<Real>(nverts * 2);
  auto ntris = nquads * 4;
  auto qv2v = mesh->ask_elem_verts();
  auto tv2v = Write<LO>(ntris * 3);
  auto f = OMEGA_H_LAMBDA(LO q) {
    auto qqv2v = gather_verts<8>(qv2v, q);
    for (Int qqt = 0; qqt < 4; ++qqt) {
      auto t = q * 4 + qqt;
      auto v0 = qqv2v[qqt];
      auto v1 = qqv2v[(qqt + 1) % 4];
      auto vq = nold_verts + q;
      tv2v[t * 3 + 0] = vq;
      tv2v[t * 3 + 1] = v0;
      tv2v[t * 3 + 2] = v1;
    }
  };
  parallel_for(nquads, f);
  map_into_range(old_coords, 0, nold_verts, coords, 2);
  map_into_range(quad_center_coords, nold_verts, nold_verts + nquads, coords, 2);
  Mesh new_mesh(mesh->library());
  build_from_elems_and_coords(&new_mesh, OMEGA_H_SIMPLEX, 2, tv2v, coords);
  *mesh = new_mesh;
}

void hexes2tets(Mesh* mesh) {
  auto nold_verts = mesh->nverts();
  auto nquads = mesh->nfaces();
  auto nhexes = mesh->nregions();
  auto old_coords = mesh->coords();
  auto quad_center_coords = average_field(mesh, 2, 3, old_coords);
  auto hex_center_coords = average_field(mesh, 3, 3, old_coords);
  auto nverts = nold_verts + nquads + nhexes;
  auto coords = Write<Real>(nverts * 3);
  auto ntets = nhexes * 6 * 4;
  auto hv2v = mesh->ask_elem_verts();
  auto hq2q = mesh->ask_down(3, 2).ab2b;
  auto tv2v = Write<LO>(ntets * 4);
  auto f = OMEGA_H_LAMBDA(LO h) {
    auto hhv2v = gather_verts<8>(hv2v, h);
    auto hhq2q = gather_down<6>(hq2q, h);
    for (Int hhq = 0; hhq < 6; ++hhq) {
      auto q = hhq2q[hhq];
      for (Int hhqt = 0; hhqt < 4; ++hhqt) {
        auto t = (h * 6 + hhq) * 4 + hhqt;
        auto v0 = hhv2v[hypercube_down_template(3, 2, hhq, hhqt)];
        auto v1 = hhv2v[hypercube_down_template(3, 2, hhq, (hhqt + 1) % 4)];
        auto vq = nold_verts + q;
        auto vh = nold_verts + nquads + h;
        tv2v[t * 4 + 0] = vq;
        tv2v[t * 4 + 1] = v1;
        tv2v[t * 4 + 2] = v0;
        tv2v[t * 4 + 3] = vh;
      }
    }
  };
  parallel_for(nhexes, f);
  map_into_range(old_coords, 0, nold_verts, coords, 3);
  map_into_range(quad_center_coords, nold_verts, nold_verts + nquads, coords, 3);
  map_into_range(hex_center_coords, nold_verts + nquads, nverts, coords, 3);
  Mesh new_mesh(mesh->library());
  build_from_elems_and_coords(&new_mesh, OMEGA_H_SIMPLEX, 3, tv2v, coords);
  *mesh = new_mesh;
}

}  // end namespace Omega_h
