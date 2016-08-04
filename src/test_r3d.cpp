#include "omega_h_r3d.hpp"
#include "algebra.hpp"

#include <iostream>

int main() {
  osh::Vector<3> verts[4] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  osh::r3d::Plane<3> faces[4];
  osh::r3d::tet_faces_from_verts(faces, verts);
  OSH_CHECK(osh::are_close(faces[0].n, -osh::normalize(osh::vector_3(1,1,1))));
  OSH_CHECK(osh::are_close(faces[0].d, -faces[0].n[0]));
  for (osh::Int i = 0; i < 3; ++i) {
    auto v = osh::vector_3(0,0,0);
    v[i] = 1;
    OSH_CHECK(osh::are_close(faces[i + 1].n, v));
    OSH_CHECK(osh::are_close(faces[i + 1].d, 0));
  }
  osh::r3d::Polytope<3> a;
  osh::r3d::init_tet(&a, verts);
  OSH_CHECK(a.nverts == 4);
  constexpr osh::Int nplanes = 1;
  osh::r3d::Plane<3> planes[nplanes] = {
    {{1,0,0},-0.5}
  };
  osh::r3d::clip(&a, planes, nplanes);
  OSH_CHECK(a.nverts == 4);
  auto volume = osh::r3d::integrate(&a,
      osh::r3d::Polynomial<3,0>{{1}});
  OSH_CHECK(osh::are_close(volume, osh::cube(0.5) / 6.0));
}
