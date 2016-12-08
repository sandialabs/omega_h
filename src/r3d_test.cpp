#include "Omega_h_r3d.hpp"
#include "algebra.hpp"

#include <iostream>

static void test_3d() {
  Omega_h::Few<Omega_h::Vector<3>, 4> verts = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  auto faces = Omega_h::r3d::faces_from_verts(verts);
  OMEGA_H_CHECK(Omega_h::are_close(
      faces[0].n, -Omega_h::normalize(Omega_h::vector_3(1, 1, 1))));
  OMEGA_H_CHECK(Omega_h::are_close(faces[0].d, -faces[0].n[0]));
  for (Omega_h::Int i = 0; i < 3; ++i) {
    auto v = Omega_h::vector_3(0, 0, 0);
    v[i] = 1;
    OMEGA_H_CHECK(Omega_h::are_close(faces[i + 1].n, v));
    OMEGA_H_CHECK(Omega_h::are_close(faces[i + 1].d, 0));
  }
  Omega_h::r3d::Polytope<3> a;
  Omega_h::r3d::init(a, verts);
  OMEGA_H_CHECK(a.nverts == 4);
  auto b = a;
  Omega_h::r3d::clip(b, Omega_h::Few<Omega_h::r3d::Plane<3>, 1>({{{1, 0, 0}, -0.5}}));
  OMEGA_H_CHECK(b.nverts == 4);
  auto volume = Omega_h::r3d::measure(b);
  OMEGA_H_CHECK(Omega_h::are_close(volume, Omega_h::cube(0.5) / 6.0));
  auto c = a;
  Omega_h::r3d::clip(c, Omega_h::Few<Omega_h::r3d::Plane<3>, 1>({{{-1, 0, 0}, 0.5}}));
  OMEGA_H_CHECK(c.nverts == 6);
  volume = Omega_h::r3d::measure(c);
  OMEGA_H_CHECK(
      Omega_h::are_close(volume, (1. / 6.) - (Omega_h::cube(0.5) / 6.)));
  Omega_h::Few<Omega_h::Vector<3>, 4> verts2 = {
      {1, 0, 0}, {1, 1, 0}, {0, 0, 0}, {1, 0, 1}};
  Omega_h::r3d::Polytope<3> intersection;
  Omega_h::r3d::intersect_simplices(intersection, verts, verts2);
  volume = Omega_h::r3d::measure(intersection);
  OMEGA_H_CHECK(Omega_h::are_close(volume, (1. / 3.) * (1. / 4.) * (1. / 2.)));
  Omega_h::Few<Omega_h::Vector<3>, 4> far_verts = {
      {0, 0, 4}, {1, 0, 4}, {0, 1, 4}, {0, 0, 5}};
  Omega_h::r3d::Polytope<3> null_intersection;
  Omega_h::r3d::intersect_simplices(null_intersection, verts, far_verts);
  OMEGA_H_CHECK(null_intersection.nverts == 0);
}

static void test_2d() {
  Omega_h::Few<Omega_h::Vector<2>, 3> verts = {{0, 0}, {1, 0}, {0, 1}};
  auto faces = Omega_h::r3d::faces_from_verts(verts);
  OMEGA_H_CHECK(Omega_h::are_close(faces[0].n, Omega_h::vector_2(0, 1)));
  OMEGA_H_CHECK(Omega_h::are_close(faces[0].d, 0));
  OMEGA_H_CHECK(Omega_h::are_close(
      faces[1].n, -Omega_h::normalize(Omega_h::vector_2(1, 1))));
  OMEGA_H_CHECK(Omega_h::are_close(faces[1].d, -faces[1].n[0]));
  OMEGA_H_CHECK(Omega_h::are_close(faces[2].n, Omega_h::vector_2(1, 0)));
  OMEGA_H_CHECK(Omega_h::are_close(faces[2].d, 0));
  Omega_h::r3d::Polytope<2> a;
  Omega_h::r3d::init(a, verts);
  OMEGA_H_CHECK(a.nverts == 3);
  auto b = a;
  Omega_h::r3d::clip(b, Omega_h::Few<Omega_h::r3d::Plane<2>, 1>({{{1, 0}, -0.5}}));
  OMEGA_H_CHECK(b.nverts == 3);
  auto area = Omega_h::r3d::measure(b);
  OMEGA_H_CHECK(Omega_h::are_close(area, Omega_h::square(0.5) / 2.0));
  auto c = a;
  Omega_h::r3d::clip(c, Omega_h::Few<Omega_h::r3d::Plane<2>, 1>({{{-1, 0}, 0.5}}));
  OMEGA_H_CHECK(c.nverts == 4);
  area = Omega_h::r3d::measure(c);
  OMEGA_H_CHECK(
      Omega_h::are_close(area, (1. / 2.) - (Omega_h::square(0.5) / 2.0)));
  Omega_h::Few<Omega_h::Vector<2>, 3> verts2 = {{0, 0}, {1, 0}, {1, 1}};
  Omega_h::r3d::Polytope<2> intersection;
  Omega_h::r3d::intersect_simplices(intersection, verts, verts2);
  area = Omega_h::r3d::measure(intersection);
  OMEGA_H_CHECK(Omega_h::are_close(area, 1. / 4.));
  Omega_h::Few<Omega_h::Vector<2>, 3> far_verts = {{0, 4}, {1, 4}, {0, 5}};
  Omega_h::r3d::Polytope<2> null_intersection;
  Omega_h::r3d::intersect_simplices(null_intersection, verts, far_verts);
  OMEGA_H_CHECK(null_intersection.nverts == 0);
}

int main() {
  test_2d();
  test_3d();
}
