#include "algebra.hpp"
#include "omega_h_r3d.hpp"

#include <iostream>

static void test_3d() {
  osh::Few<osh::Vector<3>, 4> verts = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  auto faces = osh::r3d::faces_from_verts(verts);
  OSH_CHECK(
      osh::are_close(faces[0].n, -osh::normalize(osh::vector_3(1, 1, 1))));
  OSH_CHECK(osh::are_close(faces[0].d, -faces[0].n[0]));
  for (osh::Int i = 0; i < 3; ++i) {
    auto v = osh::vector_3(0, 0, 0);
    v[i] = 1;
    OSH_CHECK(osh::are_close(faces[i + 1].n, v));
    OSH_CHECK(osh::are_close(faces[i + 1].d, 0));
  }
  auto a = osh::r3d::init(verts);
  OSH_CHECK(a.nverts == 4);
  auto b =
      osh::r3d::clip(a, osh::Few<osh::r3d::Plane<3>, 1>({{{1, 0, 0}, -0.5}}));
  OSH_CHECK(b.nverts == 4);
  auto volume = osh::r3d::measure(b);
  OSH_CHECK(osh::are_close(volume, osh::cube(0.5) / 6.0));
  auto c =
      osh::r3d::clip(a, osh::Few<osh::r3d::Plane<3>, 1>({{{-1, 0, 0}, 0.5}}));
  OSH_CHECK(c.nverts == 6);
  volume = osh::r3d::measure(c);
  OSH_CHECK(osh::are_close(volume, (1. / 6.) - (osh::cube(0.5) / 6.)));
}

static void test_2d() {
  osh::Few<osh::Vector<2>, 3> verts = {{0, 0}, {1, 0}, {0, 1}};
  auto faces = osh::r3d::faces_from_verts(verts);
  OSH_CHECK(osh::are_close(faces[0].n, osh::vector_2(0, 1)));
  OSH_CHECK(osh::are_close(faces[0].d, 0));
  OSH_CHECK(osh::are_close(faces[1].n, -osh::normalize(osh::vector_2(1, 1))));
  OSH_CHECK(osh::are_close(faces[1].d, -faces[1].n[0]));
  OSH_CHECK(osh::are_close(faces[2].n, osh::vector_2(1, 0)));
  OSH_CHECK(osh::are_close(faces[2].d, 0));
  auto a = osh::r3d::init(verts);
  OSH_CHECK(a.nverts == 3);
  auto b = osh::r3d::clip(a, osh::Few<osh::r3d::Plane<2>, 1>({{{1, 0}, -0.5}}));
  OSH_CHECK(b.nverts == 3);
  auto area = osh::r3d::measure(b);
  OSH_CHECK(osh::are_close(area, osh::square(0.5) / 2.0));
  auto c = osh::r3d::clip(a, osh::Few<osh::r3d::Plane<2>, 1>({{{-1, 0}, 0.5}}));
  OSH_CHECK(c.nverts == 4);
  area = osh::r3d::measure(c);
  OSH_CHECK(osh::are_close(area, (1. / 2.) - (osh::square(0.5) / 2.0)));
}

int main() {
  test_3d();
  test_2d();
}
