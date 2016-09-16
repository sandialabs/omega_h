#include "Omega_h_r3d.hpp"
#include "algebra.hpp"
#include "basis_polynomial.hpp"

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
  auto a = Omega_h::r3d::init(verts);
  OMEGA_H_CHECK(a.nverts == 4);
  auto b = Omega_h::r3d::clip(
      a, Omega_h::Few<Omega_h::r3d::Plane<3>, 1>({{{1, 0, 0}, -0.5}}));
  OMEGA_H_CHECK(b.nverts == 4);
  auto volume = Omega_h::r3d::measure(b);
  OMEGA_H_CHECK(Omega_h::are_close(volume, Omega_h::cube(0.5) / 6.0));
  auto c = Omega_h::r3d::clip(
      a, Omega_h::Few<Omega_h::r3d::Plane<3>, 1>({{{-1, 0, 0}, 0.5}}));
  OMEGA_H_CHECK(c.nverts == 6);
  volume = Omega_h::r3d::measure(c);
  OMEGA_H_CHECK(
      Omega_h::are_close(volume, (1. / 6.) - (Omega_h::cube(0.5) / 6.)));
  Omega_h::Few<Omega_h::Vector<3>, 4> verts2 = {
      {1, 0, 0}, {1, 1, 0}, {0, 0, 0}, {1, 0, 1}};
  volume =
      Omega_h::r3d::measure(Omega_h::r3d::intersect_simplices(verts, verts2));
  OMEGA_H_CHECK(Omega_h::are_close(volume, (1. / 3.) * (1. / 4.) * (1. / 2.)));
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
  auto a = Omega_h::r3d::init(verts);
  OMEGA_H_CHECK(a.nverts == 3);
  auto b = Omega_h::r3d::clip(
      a, Omega_h::Few<Omega_h::r3d::Plane<2>, 1>({{{1, 0}, -0.5}}));
  OMEGA_H_CHECK(b.nverts == 3);
  auto area = Omega_h::r3d::measure(b);
  OMEGA_H_CHECK(Omega_h::are_close(area, Omega_h::square(0.5) / 2.0));
  auto c = Omega_h::r3d::clip(
      a, Omega_h::Few<Omega_h::r3d::Plane<2>, 1>({{{-1, 0}, 0.5}}));
  OMEGA_H_CHECK(c.nverts == 4);
  area = Omega_h::r3d::measure(c);
  OMEGA_H_CHECK(
      Omega_h::are_close(area, (1. / 2.) - (Omega_h::square(0.5) / 2.0)));
  Omega_h::Few<Omega_h::Vector<2>, 3> verts2 = {{0, 0}, {1, 0}, {1, 1}};
  area =
      Omega_h::r3d::measure(Omega_h::r3d::intersect_simplices(verts, verts2));
  OMEGA_H_CHECK(Omega_h::are_close(area, 1. / 4.));
}

template <Omega_h::Int dim>
static void test_pair_integral_dim(
    Omega_h::Few<Omega_h::Vector<dim>, dim + 1> elem_pts) {
  std::cerr << "TEST_PAIR_INTEGRAL_DIM<" << dim << ">\n";
  auto size = Omega_h::element_size(Omega_h::simplex_basis<dim, dim>(elem_pts));
  auto polytope = Omega_h::r3d::init(elem_pts);
  for (Omega_h::Int i = 0; i <= dim; ++i) {
    auto polynomial1 = Omega_h::get_basis_polynomial(elem_pts, i);
    std::cout << "first polynomial";
    for (Omega_h::Int j = 0; j <= dim; ++j)
      std::cout << ' ' << polynomial1.coeffs[j];
    std::cout << '\n';
    for (Omega_h::Int j = 0; j <= dim; ++j) {
      if (i == j) continue;
      auto polynomial2 = Omega_h::get_basis_polynomial(elem_pts, j);
      std::cout << "second polynomial";
      for (Omega_h::Int k = 0; k <= dim; ++k)
        std::cout << ' ' << polynomial2.coeffs[k];
      std::cout << '\n';
      auto pair_polynomial = polynomial1 * polynomial2;
      auto integral = Omega_h::r3d::integrate(polytope, pair_polynomial);
      std::cout << "integral " << integral << '\n';
      OMEGA_H_CHECK(
          Omega_h::are_close(integral, size / ((dim + 1) * (dim + 2))));
    }
  }
}

static void test_pair_integrals() {
  Omega_h::Few<Omega_h::Vector<2>, 3> parent_tri = {
      {0, 0}, {1, 0}, {0, 1},
  };
  test_pair_integral_dim<2>(parent_tri);
  Omega_h::Few<Omega_h::Vector<3>, 4> parent_tet = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
  };
  test_pair_integral_dim<3>(parent_tet);
  // failing, will fix soon
  // Omega_h::Few<Omega_h::Vector<2>, 3> perfect_tri(
  //    {{1, 0}, {0, sqrt(3.0)}, {-1, 0}});
  // test_pair_integral_dim<2>(perfect_tri);
  // Omega_h::Few<Omega_h::Vector<3>, 4> perfect_tet(
  //    {{1, 0, -1.0 / sqrt(2.0)}, {1, 0, -1.0 / sqrt(2.0)},
  //     {0, -1, 1.0 / sqrt(2.0)}, {0, 1, 1.0 / sqrt(2.0)}});
  // test_pair_integral_dim<3>(perfect_tet);
}

int main() {
  test_3d();
  test_2d();
  test_pair_integrals();
}
