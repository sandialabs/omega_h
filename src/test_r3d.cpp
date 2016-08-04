#include "omega_h_r3d.hpp"
#include "algebra.hpp"

#include <iostream>

int main() {
  osh::Few<osh::Vector<3>, 4> verts = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  auto faces = osh::r3d::tet_faces_from_verts(verts);
  OSH_CHECK(osh::are_close(faces[0].n, -osh::normalize(osh::vector_3(1,1,1))));
  OSH_CHECK(osh::are_close(faces[0].d, -faces[0].n[0]));
  for (osh::Int i = 0; i < 3; ++i) {
    auto v = osh::vector_3(0,0,0);
    v[i] = 1;
    OSH_CHECK(osh::are_close(faces[i + 1].n, v));
    OSH_CHECK(osh::are_close(faces[i + 1].d, 0));
  }
  auto a = osh::r3d::init_tet(verts);
  OSH_CHECK(a.nverts == 4);
  auto b = osh::r3d::clip(a, osh::Few<osh::r3d::Plane<3>, 1>({{{1,0,0},-0.5}}));
  OSH_CHECK(b.nverts == 4);
  auto volume = osh::r3d::measure(b);
  OSH_CHECK(osh::are_close(volume, osh::cube(0.5) / 6.0));
  auto c = osh::r3d::clip(a, osh::Few<osh::r3d::Plane<3>, 1>({{{-1,0,0},0.5}}));
  OSH_CHECK(c.nverts == 6);
  volume = osh::r3d::measure(c);
  OSH_CHECK(osh::are_close(volume, (1. / 6.) - (osh::cube(0.5) / 6.)));
}
