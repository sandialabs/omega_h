#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_bbox.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_confined.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_hilbert.hpp"
#include "Omega_h_hypercube.hpp"
#include "Omega_h_inertia.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_recover.hpp"
#include "Omega_h_refine_qualities.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_swap2d.hpp"
#include "Omega_h_swap3d_choice.hpp"
#include "Omega_h_swap3d_loop.hpp"

#include <sstream>

using namespace Omega_h;

static void test_down_template() {
  OMEGA_H_CHECK(simplex_down_template(1, 0, 0, 0) == 0);
  OMEGA_H_CHECK(simplex_down_template(1, 0, 1, 0) == 1);
  OMEGA_H_CHECK(simplex_down_template(2, 0, 0, 0) == 0);
  OMEGA_H_CHECK(simplex_down_template(2, 0, 1, 0) == 1);
  OMEGA_H_CHECK(simplex_down_template(2, 0, 2, 0) == 2);
  OMEGA_H_CHECK(simplex_down_template(2, 1, 0, 0) == 0);
  OMEGA_H_CHECK(simplex_down_template(2, 1, 0, 1) == 1);
  OMEGA_H_CHECK(simplex_down_template(2, 1, 1, 0) == 1);
  OMEGA_H_CHECK(simplex_down_template(2, 1, 1, 1) == 2);
  OMEGA_H_CHECK(simplex_down_template(2, 1, 2, 0) == 2);
  OMEGA_H_CHECK(simplex_down_template(2, 1, 2, 1) == 0);
  OMEGA_H_CHECK(simplex_down_template(3, 0, 0, 0) == 0);
  OMEGA_H_CHECK(simplex_down_template(3, 0, 1, 0) == 1);
  OMEGA_H_CHECK(simplex_down_template(3, 0, 2, 0) == 2);
  OMEGA_H_CHECK(simplex_down_template(3, 0, 3, 0) == 3);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 0, 0) == 0);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 0, 1) == 1);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 1, 0) == 1);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 1, 1) == 2);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 2, 0) == 2);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 2, 1) == 0);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 3, 0) == 0);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 3, 1) == 3);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 4, 0) == 1);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 4, 1) == 3);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 5, 0) == 2);
  OMEGA_H_CHECK(simplex_down_template(3, 1, 5, 1) == 3);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 0, 0) == 0);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 0, 1) == 2);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 0, 2) == 1);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 1, 0) == 0);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 1, 1) == 1);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 1, 2) == 3);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 2, 0) == 1);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 2, 1) == 2);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 2, 2) == 3);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 3, 0) == 2);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 3, 1) == 0);
  OMEGA_H_CHECK(simplex_down_template(3, 2, 3, 2) == 3);
}

static OMEGA_H_DEVICE bool same_adj(Int a[], Int b[]) {
  for (Int i = 0; i < 3; ++i)
    if (a[i] != b[i]) return false;
  return true;
}

void test_tri_align() {
  auto f = OMEGA_H_LAMBDA(LO) {
    Int ident[3] = {0, 1, 2};
    Int out[3];
    Int out2[3];
    /* check that flipping and rotating do what we want */
    {
      align_adj<3>(make_code(true, 0, 0), ident, 0, out, 0);
      Int expect[3] = {0, 2, 1};
      OMEGA_H_CHECK(same_adj(out, expect));
    }
    {
      align_adj<3>(make_code(false, 1, 0), ident, 0, out, 0);
      Int expect[3] = {2, 0, 1};
      OMEGA_H_CHECK(same_adj(out, expect));
    }
    {
      align_adj<3>(make_code(false, 2, 0), ident, 0, out, 0);
      Int expect[3] = {1, 2, 0};
      OMEGA_H_CHECK(same_adj(out, expect));
    }
    /* check that compound_alignments does its job */
    for (I8 rot1 = 0; rot1 < 3; ++rot1)
      for (I8 flip1 = 0; flip1 < 2; ++flip1)
        for (I8 rot2 = 0; rot2 < 3; ++rot2)
          for (I8 flip2 = 0; flip2 < 2; ++flip2) {
            I8 code1 = make_code(flip1, rot1, 0);
            I8 code2 = make_code(flip2, rot2, 0);
            align_adj<3>(code1, ident, 0, out, 0);
            align_adj<3>(code2, out, 0, out2, 0);
            Int out3[3];
            I8 code3 = compound_alignments<3>(code1, code2);
            align_adj<3>(code3, ident, 0, out3, 0);
            OMEGA_H_CHECK(same_adj(out2, out3));
          }
  };
  parallel_for(1, f);
}

static void test_form_uses() {
  OMEGA_H_CHECK(form_uses(LOs({0, 1, 2}), OMEGA_H_SIMPLEX, 2, 1) ==
                LOs({0, 1, 1, 2, 2, 0}));
  OMEGA_H_CHECK(form_uses(LOs({0, 1, 2, 3}), OMEGA_H_SIMPLEX, 3, 1) ==
                LOs({0, 1, 1, 2, 2, 0, 0, 3, 1, 3, 2, 3}));
  OMEGA_H_CHECK(form_uses(LOs({0, 1, 2, 3}), OMEGA_H_SIMPLEX, 3, 2) ==
                LOs({0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3}));
  OMEGA_H_CHECK(form_uses(LOs({0, 1, 2, 3}), OMEGA_H_HYPERCUBE, 2, 1) ==
                LOs({0, 1, 1, 2, 2, 3, 3, 0}));
  OMEGA_H_CHECK(form_uses(LOs({0, 1, 2, 3, 4, 5, 6, 7}), OMEGA_H_HYPERCUBE, 3,
                    1) == LOs({0, 1, 1, 2, 2, 3, 3, 0, 0, 4, 1, 5, 2, 6, 3, 7,
                              4, 5, 5, 6, 6, 7, 7, 4}));
  OMEGA_H_CHECK(form_uses(LOs({0, 1, 2, 3, 4, 5, 6, 7}), OMEGA_H_HYPERCUBE, 3,
                    2) == LOs({1, 0, 3, 2, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6,
                              3, 0, 4, 7, 4, 5, 6, 7}));
}

static void test_reflect_down() {
  Adj a;
  a = reflect_down(LOs({}), LOs({}), OMEGA_H_SIMPLEX, 0, 2, 1);
  OMEGA_H_CHECK(a.ab2b == LOs({}));
  OMEGA_H_CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({}), LOs({}), OMEGA_H_SIMPLEX, 0, 3, 1);
  OMEGA_H_CHECK(a.ab2b == LOs({}));
  OMEGA_H_CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({}), LOs({}), OMEGA_H_SIMPLEX, 0, 3, 2);
  OMEGA_H_CHECK(a.ab2b == LOs({}));
  OMEGA_H_CHECK(a.codes == Read<I8>({}));
  a = reflect_down(
      LOs({0, 1, 2}), LOs({0, 1, 1, 2, 2, 0}), OMEGA_H_SIMPLEX, 3, 2, 1);
  OMEGA_H_CHECK(a.ab2b == LOs({0, 1, 2}));
  OMEGA_H_CHECK(a.codes == Read<I8>({0, 0, 0}));
  a = reflect_down(LOs({0, 1, 2, 3}), LOs({0, 1, 1, 2, 2, 0, 0, 3, 1, 3, 2, 3}),
      OMEGA_H_SIMPLEX, 4, 3, 1);
  OMEGA_H_CHECK(a.ab2b == LOs({0, 1, 2, 3, 4, 5}));
  OMEGA_H_CHECK(a.codes == Read<I8>({0, 0, 0, 0, 0, 0}));
  a = reflect_down(LOs({0, 1, 2, 3}), LOs({0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3}),
      OMEGA_H_SIMPLEX, 4, 3, 2);
  OMEGA_H_CHECK(a.ab2b == LOs({0, 1, 2, 3}));
  OMEGA_H_CHECK(a.codes == Read<I8>({0, 0, 0, 0}));
  a = reflect_down(LOs({0, 1, 2, 3}), LOs({0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0}),
      OMEGA_H_SIMPLEX, 4, 3, 2);
  OMEGA_H_CHECK(a.ab2b == LOs({0, 1, 2, 3}));
  OMEGA_H_CHECK(a.codes == Read<I8>(4, make_code(true, 0, 0)));
  a = reflect_down(LOs({0, 1, 2, 2, 3, 0}), LOs({0, 1, 1, 2, 2, 3, 3, 0, 0, 2}),
      OMEGA_H_SIMPLEX, 4, 2, 1);
  OMEGA_H_CHECK(a.ab2b == LOs({0, 1, 4, 2, 3, 4}));
  a = reflect_down(LOs({}), LOs({}), OMEGA_H_HYPERCUBE, 0, 2, 1);
  OMEGA_H_CHECK(a.ab2b == LOs({}));
  OMEGA_H_CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({}), LOs({}), OMEGA_H_HYPERCUBE, 0, 3, 1);
  OMEGA_H_CHECK(a.ab2b == LOs({}));
  OMEGA_H_CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({}), LOs({}), OMEGA_H_HYPERCUBE, 0, 3, 2);
  OMEGA_H_CHECK(a.ab2b == LOs({}));
  OMEGA_H_CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({0, 1, 2, 3}), LOs({0, 1, 1, 2, 2, 3, 3, 0}),
      OMEGA_H_HYPERCUBE, 4, 2, 1);
  OMEGA_H_CHECK(a.ab2b == LOs({0, 1, 2, 3}));
  OMEGA_H_CHECK(a.codes == Read<I8>({0, 0, 0, 0}));
  auto hex_verts = LOs(8, 0, 1);
  a = reflect_down(hex_verts,
      LOs({0, 1, 1, 2, 2, 3, 3, 0, 0, 4, 1, 5, 2, 6, 3, 7, 4, 5, 5, 6, 6, 7, 7,
          4}),
      OMEGA_H_HYPERCUBE, 8, 3, 1);
  OMEGA_H_CHECK(a.ab2b == LOs(12, 0, 1));
  OMEGA_H_CHECK(a.codes == Read<I8>(12, 0));
  a = reflect_down(hex_verts,
      LOs({1, 0, 3, 2, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7, 4, 5, 6,
          7}),
      OMEGA_H_HYPERCUBE, 8, 3, 2);
  OMEGA_H_CHECK(a.ab2b == LOs(6, 0, 1));
  OMEGA_H_CHECK(a.codes == Read<I8>(6, 0));
  a = reflect_down(hex_verts,
      LOs({2, 3, 0, 1, 4, 5, 1, 0, 5, 6, 2, 1, 6, 7, 3, 2, 7, 4, 0, 3, 7, 6, 5,
          4}),
      OMEGA_H_HYPERCUBE, 8, 3, 2);
  OMEGA_H_CHECK(a.ab2b == LOs(6, 0, 1));
  OMEGA_H_CHECK(a.codes == Read<I8>(6, make_code(true, 1, 0)));
  // a = reflect_down(
  //    LOs({0, 1, 2, 2, 3, 0}), LOs({0, 1, 1, 2, 2, 3, 3, 0, 0, 2}),
  //    OMEGA_H_HYPERCUBE, 4, 2, 1);
  // OMEGA_H_CHECK(a.ab2b == LOs({0, 1, 4, 2, 3, 4}));
}

static void test_find_unique() {
  OMEGA_H_CHECK(find_unique(LOs({}), OMEGA_H_SIMPLEX, 2, 1) == LOs({}));
  OMEGA_H_CHECK(find_unique(LOs({}), OMEGA_H_SIMPLEX, 3, 1) == LOs({}));
  OMEGA_H_CHECK(find_unique(LOs({}), OMEGA_H_SIMPLEX, 3, 2) == LOs({}));
  OMEGA_H_CHECK(find_unique(LOs({0, 1, 2, 2, 3, 0}), OMEGA_H_SIMPLEX, 2, 1) ==
                LOs({0, 1, 0, 2, 3, 0, 1, 2, 2, 3}));
  OMEGA_H_CHECK(find_unique(LOs({}), OMEGA_H_HYPERCUBE, 2, 1) == LOs({}));
  OMEGA_H_CHECK(find_unique(LOs({}), OMEGA_H_HYPERCUBE, 3, 1) == LOs({}));
  OMEGA_H_CHECK(find_unique(LOs({}), OMEGA_H_HYPERCUBE, 3, 2) == LOs({}));
  auto a = find_unique(LOs({0, 1, 2, 3}), OMEGA_H_HYPERCUBE, 2, 1);
  OMEGA_H_CHECK(find_unique(LOs({0, 1, 2, 3}), OMEGA_H_HYPERCUBE, 2, 1) ==
                LOs({0, 1, 3, 0, 1, 2, 2, 3}));
}

static void test_hilbert() {
  /* this is the original test from Skilling's paper */
  hilbert::coord_t X[3] = {5, 10, 20};  // any position in 32x32x32 cube
  hilbert::AxestoTranspose(X, 5,
      3);  // Hilbert transpose for 5 bits and 3 dimensions
  std::stringstream stream;
  stream << "Hilbert integer = " << (X[0] >> 4 & 1) << (X[1] >> 4 & 1)
         << (X[2] >> 4 & 1) << (X[0] >> 3 & 1) << (X[1] >> 3 & 1)
         << (X[2] >> 3 & 1) << (X[0] >> 2 & 1) << (X[1] >> 2 & 1)
         << (X[2] >> 2 & 1) << (X[0] >> 1 & 1) << (X[1] >> 1 & 1)
         << (X[2] >> 1 & 1) << (X[0] >> 0 & 1) << (X[1] >> 0 & 1)
         << (X[2] >> 0 & 1) << " = 7865 check";
  std::string expected = "Hilbert integer = 001111010111001 = 7865 check";
  OMEGA_H_CHECK(stream.str() == expected);
  hilbert::coord_t Y[3];
  hilbert::untranspose(X, Y, 5, 3);
  std::stringstream stream2;
  stream2 << "Hilbert integer = " << (Y[0] >> 4 & 1) << (Y[0] >> 3 & 1)
          << (Y[0] >> 2 & 1) << (Y[0] >> 1 & 1) << (Y[0] >> 0 & 1)
          << (Y[1] >> 4 & 1) << (Y[1] >> 3 & 1) << (Y[1] >> 2 & 1)
          << (Y[1] >> 1 & 1) << (Y[1] >> 0 & 1) << (Y[2] >> 4 & 1)
          << (Y[2] >> 3 & 1) << (Y[2] >> 2 & 1) << (Y[2] >> 1 & 1)
          << (Y[2] >> 0 & 1) << " = 7865 check";
  OMEGA_H_CHECK(stream2.str() == expected);
}

static void test_bbox() {
  OMEGA_H_CHECK(are_close(BBox<2>(vector_2(-3, -3), vector_2(3, 3)),
      find_bounding_box<2>(Reals({0, -3, 3, 0, 0, 3, -3, 0}))));
  OMEGA_H_CHECK(are_close(BBox<3>(vector_3(-3, -3, -3), vector_3(3, 3, 3)),
      find_bounding_box<3>(
          Reals({0, -3, 0, 3, 0, 0, 0, 3, 0, -3, 0, 0, 0, 0, -3, 0, 0, 3}))));
}

static void test_build(Library* lib) {
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 2, LOs({0, 1, 2}), 3);
    OMEGA_H_CHECK(mesh.ask_down(2, 0).ab2b == LOs({0, 1, 2}));
    OMEGA_H_CHECK(mesh.ask_down(2, 1).ab2b == LOs({0, 2, 1}));
    OMEGA_H_CHECK(mesh.ask_down(1, 0).ab2b == LOs({0, 1, 2, 0, 1, 2}));
  }
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 3, LOs({0, 1, 2, 3}), 4);
    OMEGA_H_CHECK(mesh.ask_down(3, 0).ab2b == LOs({0, 1, 2, 3}));
  }
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_HYPERCUBE, 2, LOs({0, 1, 2, 3}), 4);
    OMEGA_H_CHECK(mesh.ask_down(2, 0).ab2b == LOs({0, 1, 2, 3}));
    OMEGA_H_CHECK(mesh.ask_down(1, 0).ab2b == LOs({0, 1, 3, 0, 1, 2, 2, 3}));
    OMEGA_H_CHECK(mesh.ask_down(2, 1).ab2b == LOs({0, 2, 3, 1}));
  }
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_HYPERCUBE, 3, LOs(8, 0, 1), 8);
    OMEGA_H_CHECK(mesh.ask_down(3, 0).ab2b == LOs(8, 0, 1));
    OMEGA_H_CHECK(
        mesh.ask_down(2, 0).ab2b == LOs({1, 0, 3, 2, 0, 1, 5, 4, 3, 0, 4, 7, 1,
                                        2, 6, 5, 2, 3, 7, 6, 4, 5, 6, 7}));
    OMEGA_H_CHECK(
        mesh.ask_up(0, 2).ab2b == LOs({0, 1, 2, 0, 1, 3, 0, 3, 4, 0, 2, 4, 1, 2,
                                      5, 1, 3, 5, 3, 4, 5, 2, 4, 5}));
  }
  {
    auto mesh =
        build_box(lib->world(), OMEGA_H_HYPERCUBE, 1.0, 1.0, 1.0, 2, 2, 2);
  }
}

static void test_star(Library* lib) {
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 2, LOs({0, 1, 2}), 3);
    Adj v2v = mesh.ask_star(VERT);
    OMEGA_H_CHECK(v2v.a2ab == LOs(4, 0, 2));
    OMEGA_H_CHECK(v2v.ab2b == LOs({1, 2, 0, 2, 0, 1}));
    Adj e2e = mesh.ask_star(EDGE);
    OMEGA_H_CHECK(e2e.a2ab == LOs(4, 0, 2));
    OMEGA_H_CHECK(e2e.ab2b == LOs({2, 1, 0, 2, 1, 0}));
  }
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 3, LOs({0, 1, 2, 3}), 4);
    Adj v2v = mesh.ask_star(VERT);
    OMEGA_H_CHECK(v2v.a2ab == LOs(5, 0, 3));
    OMEGA_H_CHECK(v2v.ab2b == LOs({1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2}));
    Adj e2e = mesh.ask_star(EDGE);
    OMEGA_H_CHECK(e2e.a2ab == LOs(7, 0, 5));
    OMEGA_H_CHECK(
        e2e.ab2b == LOs({1, 3, 4, 2, 5, 3, 0, 2, 5, 4, 0, 4, 5, 1, 3, 0, 1, 5,
                        4, 2, 2, 0, 3, 5, 1, 1, 2, 4, 3, 0}));
  }
}

static void test_dual(Library* lib) {
  Mesh mesh(lib);
  build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 2, LOs({0, 1, 2, 2, 3, 0}), 4);
  auto t2t = mesh.ask_dual();
  auto t2tt = t2t.a2ab;
  auto tt2t = t2t.ab2b;
  OMEGA_H_CHECK(t2tt == offset_scan(LOs({1, 1})));
  OMEGA_H_CHECK(tt2t == LOs({1, 0}));
}

static void test_quality() {
  Few<Vector<2>, 3> perfect_tri(
      {vector_2(1, 0), vector_2(0, std::sqrt(3.0)), vector_2(-1, 0)});
  Few<Vector<3>, 4> perfect_tet({vector_3(1, 0, -1.0 / std::sqrt(2.0)),
      vector_3(-1, 0, -1.0 / std::sqrt(2.0)),
      vector_3(0, -1, 1.0 / std::sqrt(2.0)),
      vector_3(0, 1, 1.0 / std::sqrt(2.0))});
  Few<Vector<2>, 3> flat_tri({vector_2(1, 0), vector_2(0, 0), vector_2(-1, 0)});
  Few<Vector<3>, 4> flat_tet({vector_3(1, 0, 0), vector_3(-1, 0, 0),
      vector_3(0, -1, 0), vector_3(0, 1, 0)});
  Few<Vector<2>, 3> inv_tri(
      {vector_2(1, 0), vector_2(-1, 0), vector_2(0, std::sqrt(3.0))});
  Few<Vector<3>, 4> inv_tet({vector_3(1, 0, -1.0 / std::sqrt(2.0)),
      vector_3(-1, 0, -1.0 / std::sqrt(2.0)),
      vector_3(0, 1, 1.0 / std::sqrt(2.0)),
      vector_3(0, -1, 1.0 / std::sqrt(2.0))});
  Tensor<2> id_metric_2 = identity_tensor<2>();
  Tensor<3> id_metric_3 = identity_tensor<3>();
  Tensor<2> x_metric_2 =
      compose_metric(identity_matrix<2, 2>(), vector_2(1, 0.5));
  Tensor<3> x_metric_3 =
      compose_metric(identity_matrix<3, 3>(), vector_3(1, 1, 0.5));
  Few<Vector<2>, 3> x_tri;
  for (Int i = 0; i < 3; ++i) {
    x_tri[i] = perfect_tri[i];
    x_tri[i][1] /= 2;
  }
  Few<Vector<3>, 4> x_tet;
  for (Int i = 0; i < 4; ++i) {
    x_tet[i] = perfect_tet[i];
    x_tet[i][2] /= 2;
  }
  OMEGA_H_CHECK(
      are_close(metric_element_quality(perfect_tri, id_metric_2), 1.0));
  OMEGA_H_CHECK(
      are_close(metric_element_quality(perfect_tet, id_metric_3), 1.0));
  OMEGA_H_CHECK(are_close(metric_element_quality(flat_tri, id_metric_2), 0.0));
  OMEGA_H_CHECK(are_close(metric_element_quality(flat_tet, id_metric_3), 0.0));
  OMEGA_H_CHECK(metric_element_quality(inv_tri, id_metric_2) < 0.0);
  OMEGA_H_CHECK(metric_element_quality(inv_tet, id_metric_3) < 0.0);
  OMEGA_H_CHECK(are_close(metric_element_quality(x_tri, x_metric_2), 1.0));
  OMEGA_H_CHECK(are_close(metric_element_quality(x_tet, x_metric_3), 1.0));
}

static void test_inertial_bisect(Library* lib) {
  Reals coords({2, 1, 0, 2, -1, 0, -2, 1, 0, -2, -1, 0});
  Reals masses(4, 1);
  auto self = lib->self();
  Real tolerance = 0.0;
  Vector<3> axis;
  auto marked = inertia::mark_bisection(self, coords, masses, tolerance, axis);
  OMEGA_H_CHECK(marked == Read<I8>({1, 1, 0, 0}));
  marked = inertia::mark_bisection_given_axis(
      self, coords, masses, tolerance, vector_3(0, 1, 0));
  OMEGA_H_CHECK(marked == Read<I8>({1, 0, 1, 0}));
}

static void test_average_field(Library* lib) {
  auto mesh = Mesh(lib);
  build_box_internal(&mesh, OMEGA_H_SIMPLEX, 1, 1, 0, 1, 1, 0);
  Reals v2x({2, 1, 3, 2});
  auto e2x = average_field(&mesh, 2, LOs({0, 1}), 1, v2x);
  OMEGA_H_CHECK(are_close(e2x, Reals({5. / 3., 7. / 3.})));
}

static void test_refine_qualities(Library* lib) {
  auto mesh = Mesh(lib);
  build_box_internal(&mesh, OMEGA_H_SIMPLEX, 1., 1., 0., 1, 1, 0);
  LOs candidates(mesh.nedges(), 0, 1);
  mesh.add_tag(VERT, "metric", symm_ncomps(2),
      repeat_symm(mesh.nverts(), identity_matrix<2, 2>()));
  auto quals = refine_qualities(&mesh, candidates);
  OMEGA_H_CHECK(are_close(
      quals, Reals({0.494872, 0.494872, 0.866025, 0.494872, 0.494872}), 1e-4));
}

static void test_mark_up_down(Library* lib) {
  auto mesh = Mesh(lib);
  build_box_internal(&mesh, OMEGA_H_SIMPLEX, 1., 1., 0., 1, 1, 0);
  OMEGA_H_CHECK(
      mark_down(&mesh, FACE, VERT, Read<I8>({1, 0})) == Read<I8>({1, 1, 0, 1}));
  OMEGA_H_CHECK(
      mark_up(&mesh, VERT, FACE, Read<I8>({0, 1, 0, 0})) == Read<I8>({1, 0}));
}

static void test_compare_meshes(Library* lib) {
  auto a = build_box(lib->world(), OMEGA_H_SIMPLEX, 1., 1., 0., 4, 4, 0);
  OMEGA_H_CHECK(a == a);
  Mesh b = a;
  OMEGA_H_CHECK(a == b);
  b.add_tag<I8>(VERT, "foo", 1, Read<I8>(b.nverts(), 1));
  OMEGA_H_CHECK(!(a == b));
}

static void test_swap2d_topology(Library* lib) {
  auto mesh = Mesh(lib);
  build_box_internal(&mesh, OMEGA_H_SIMPLEX, 1., 1., 0., 1, 1, 0);
  HostFew<LOs, 3> keys2prods;
  HostFew<LOs, 3> prod_verts2verts;
  auto keys2edges = LOs({2});
  swap2d_topology(&mesh, keys2edges, &keys2prods, &prod_verts2verts);
  OMEGA_H_CHECK(prod_verts2verts[EDGE] == LOs({2, 1}));
  OMEGA_H_CHECK(prod_verts2verts[FACE] == LOs({3, 2, 1, 0, 1, 2}));
  OMEGA_H_CHECK(keys2prods[EDGE] == offset_scan(LOs({1})));
  OMEGA_H_CHECK(keys2prods[FACE] == offset_scan(LOs({2})));
}

void test_swap3d_loop(Library* lib) {
  auto mesh = Mesh(lib);
  build_box_internal(&mesh, OMEGA_H_SIMPLEX, 1, 1, 1, 1, 1, 1);
  auto edges2tets = mesh.ask_up(EDGE, REGION);
  auto edges2edge_tets = edges2tets.a2ab;
  auto edge_tets2tets = edges2tets.ab2b;
  auto edge_tet_codes = edges2tets.codes;
  auto edge_verts2verts = mesh.ask_verts_of(EDGE);
  auto tet_verts2verts = mesh.ask_verts_of(REGION);
  auto f = OMEGA_H_LAMBDA(LO foo) {
    (void)foo;
    LO edge = 6;
    auto loop = swap3d::find_loop(edges2edge_tets, edge_tets2tets,
        edge_tet_codes, edge_verts2verts, tet_verts2verts, edge);
    OMEGA_H_CHECK(loop.eev2v[0] == 7);
    OMEGA_H_CHECK(loop.eev2v[1] == 0);
    OMEGA_H_CHECK(loop.size == 6);
    LO const expect[6] = {2, 3, 1, 5, 4, 6};
    for (Int i = 0; i < 6; ++i) {
      OMEGA_H_CHECK(loop.loop_verts2verts[i] == expect[i]);
    }
  };
  parallel_for(LO(1), f);
}

static void test_element_implied_metric() {
  /* perfect tri with edge lengths = 2 */
  Few<Vector<2>, 3> perfect_tri(
      {vector_2(1, 0), vector_2(0, std::sqrt(3.0)), vector_2(-1, 0)});
  auto afm = element_implied_metric(perfect_tri);
  auto bfm = compose_metric(identity_matrix<2, 2>(), vector_2(2, 2));
  OMEGA_H_CHECK(are_close(afm, bfm));
  /* perfect tet with edge lengths = 2 */
  Few<Vector<3>, 4> perfect_tet({vector_3(1, 0, -1.0 / std::sqrt(2.0)),
      vector_3(-1, 0, -1.0 / std::sqrt(2.0)),
      vector_3(0, -1, 1.0 / std::sqrt(2.0)),
      vector_3(0, 1, 1.0 / std::sqrt(2.0))});
  auto arm = element_implied_metric(perfect_tet);
  auto brm = compose_metric(identity_matrix<3, 3>(), vector_3(2, 2, 2));
  OMEGA_H_CHECK(are_close(arm, brm));
}

template <Int dim>
void test_recover_hessians_dim(Library* lib) {
  auto one_if_3d = ((dim == 3) ? 1 : 0);
  auto mesh = build_box(
      lib->world(), OMEGA_H_SIMPLEX, 1., 1., one_if_3d, 4, 4, 4 * one_if_3d);
  auto u_w = Write<Real>(mesh.nverts());
  auto coords = mesh.coords();
  // attach a field = x^2 + y^2 (+ z^2)
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto x = get_vector<dim>(coords, v);
    u_w[v] = norm_squared(x);
  };
  parallel_for(mesh.nverts(), f);
  auto u = Omega_h::Reals(u_w);
  mesh.add_tag(Omega_h::VERT, "u", 1, u);
  auto hess = recover_hessians(&mesh, u);
  // its second derivative is exactly 2dx + 2dy,
  // and both recovery steps are linear so the current
  // algorithm should get an exact answer
  Vector<dim> dv;
  for (Int i = 0; i < dim; ++i) dv[i] = 2;
  auto expected_hess = repeat_symm(mesh.nverts(), diagonal(dv));
  OMEGA_H_CHECK(are_close(hess, expected_hess));
}

static void test_recover_hessians(Library* lib) {
  test_recover_hessians_dim<2>(lib);
  test_recover_hessians_dim<3>(lib);
}

template <Int dim>
static void test_sf_scale_dim(Library* lib) {
  auto nl = 2;
  Int one_if_2d = ((dim >= 2) ? 1 : 0);
  Int one_if_3d = ((dim >= 3) ? 1 : 0);
  auto mesh = build_box(lib->world(), OMEGA_H_SIMPLEX, 1, one_if_2d, one_if_3d,
      nl, nl * one_if_2d, nl * one_if_3d);
  auto target_nelems = mesh.nelems();
  auto metrics = Omega_h::get_implied_metrics(&mesh);
  {
    // auto isos = apply_isotropy(mesh.nverts(), metrics, OMEGA_H_ISO_SIZE);
    auto isos = Omega_h::get_implied_isos(&mesh);
    auto iso_scal = get_metric_scalar_for_nelems(&mesh, isos, target_nelems);
    OMEGA_H_CHECK(are_close(iso_scal, 1.0));
  }
  {
    auto aniso_scal =
        get_metric_scalar_for_nelems(&mesh, metrics, target_nelems);
    OMEGA_H_CHECK(are_close(aniso_scal, 1.0));
  }
}

static void test_sf_scale(Library* lib) {
  test_sf_scale_dim<1>(lib);
  test_sf_scale_dim<2>(lib);
  test_sf_scale_dim<3>(lib);
}

static void test_proximity(Library* lib) {
  {  // triangle with one bridge
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 2, LOs({0, 1, 2}), 3);
    mesh.add_tag(VERT, "coordinates", 2, Reals({0, 0, 1, 0, 0, 1}));
    auto dists = get_pad_dists(&mesh, 2, Read<I8>({0, 1, 0}));
    OMEGA_H_CHECK(dists == Reals({-1.0}));
  }
  {  // triangle off-center
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 2, LOs({0, 1, 2}), 3);
    mesh.add_tag(VERT, "coordinates", 2, Reals({0, 0, 1, 1, 1, 2}));
    auto dists = get_pad_dists(&mesh, 2, Read<I8>({1, 1, 0}));
    OMEGA_H_CHECK(dists == Reals({-1.0}));
  }
  {  // triangle expected
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 2, LOs({0, 1, 2}), 3);
    mesh.add_tag(VERT, "coordinates", 2, Reals({0, 0, 1, -1, 1, 1}));
    auto dists = get_pad_dists(&mesh, 2, Read<I8>({1, 1, 0}));
    OMEGA_H_CHECK(are_close(dists, Reals({1.0})));
  }
  {  // tet with two bridges
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(VERT, "coordinates", 3, Reals(3 * 4, 0.0));
    auto dists = get_pad_dists(&mesh, 3, Read<I8>({1, 1, 0, 0, 0, 0}));
    OMEGA_H_CHECK(are_close(dists, Reals({-1.0})));
  }
  {  // tet with three bridges, off-center
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(
        VERT, "coordinates", 3, Reals({0, 0, 0, 1, -1, 1, 1, 1, 1, 1, 0, 2}));
    auto dists = get_pad_dists(&mesh, 3, Read<I8>({1, 1, 1, 0, 0, 0}));
    OMEGA_H_CHECK(are_close(dists, Reals({-1.0})));
  }
  {  // tet with three bridges, expected
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(
        VERT, "coordinates", 3, Reals({0, 0, 0, 1, -1, -1, 1, 1, -1, 1, 0, 2}));
    auto dists = get_pad_dists(&mesh, 3, Read<I8>({1, 1, 1, 0, 0, 0}));
    OMEGA_H_CHECK(are_close(dists, Reals({1.0})));
  }
  {  // edge-edge tet, off center
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(
        VERT, "coordinates", 3, Reals({0, 0, 0, 1, 0, 0, -1, 1, 0, -1, 1, 1}));
    auto dists = get_pad_dists(&mesh, 3, Read<I8>({0, 1, 1, 1, 1, 0}));
    OMEGA_H_CHECK(are_close(dists, Reals({-1.0})));
  }
  {  // edge-edge tet, expected
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, OMEGA_H_SIMPLEX, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(
        VERT, "coordinates", 3, Reals({0, 0, 0, 2, 0, 0, 1, 1, -1, 1, 1, 1}));
    auto dists = get_pad_dists(&mesh, 3, Read<I8>({0, 1, 1, 1, 1, 0}));
    OMEGA_H_CHECK(are_close(dists, Reals({1.0})));
  }
}

static void test_1d_box(Library* lib) {
  auto mesh = build_box(lib->world(), OMEGA_H_SIMPLEX, 1, 0, 0, 4, 0, 0);
}

static bool compare_hst(Int pd, Int cd, Int wc, Int wcv, SplitVertex truth) {
  auto split_vtx = hypercube_split_template(pd, cd, wc, wcv);
  if (split_vtx.dim != truth.dim) return false;
  if (split_vtx.which_down != truth.which_down) return false;
  return true;
}

static void test_hypercube_split_template() {
  OMEGA_H_CHECK(compare_hst(1, 0, 0, 0, {1, 0}));
  OMEGA_H_CHECK(compare_hst(1, 1, 0, 0, {0, 0}));
  OMEGA_H_CHECK(compare_hst(1, 1, 0, 1, {1, 0}));
  OMEGA_H_CHECK(compare_hst(1, 1, 1, 0, {1, 0}));
  OMEGA_H_CHECK(compare_hst(1, 1, 1, 1, {0, 1}));
  OMEGA_H_CHECK(compare_hst(2, 0, 0, 0, {2, 0}));
  for (Int i = 0; i <= 4; ++i) {
    OMEGA_H_CHECK(compare_hst(2, 1, i, 0, {1, i}));
    OMEGA_H_CHECK(compare_hst(2, 1, i, 1, {2, 0}));
  }
  OMEGA_H_CHECK(compare_hst(2, 2, 0, 0, {0, 0}));
  OMEGA_H_CHECK(compare_hst(2, 2, 0, 1, {1, 0}));
  OMEGA_H_CHECK(compare_hst(2, 2, 0, 2, {2, 0}));
  OMEGA_H_CHECK(compare_hst(2, 2, 0, 3, {1, 3}));
  OMEGA_H_CHECK(compare_hst(2, 2, 1, 0, {1, 0}));
  OMEGA_H_CHECK(compare_hst(2, 2, 1, 1, {0, 1}));
  OMEGA_H_CHECK(compare_hst(2, 2, 1, 2, {1, 1}));
  OMEGA_H_CHECK(compare_hst(2, 2, 1, 3, {2, 0}));
  OMEGA_H_CHECK(compare_hst(2, 2, 2, 0, {2, 0}));
  OMEGA_H_CHECK(compare_hst(2, 2, 2, 1, {1, 1}));
  OMEGA_H_CHECK(compare_hst(2, 2, 2, 2, {0, 2}));
  OMEGA_H_CHECK(compare_hst(2, 2, 2, 3, {1, 2}));
  OMEGA_H_CHECK(compare_hst(2, 2, 3, 0, {1, 3}));
  OMEGA_H_CHECK(compare_hst(2, 2, 3, 1, {2, 0}));
  OMEGA_H_CHECK(compare_hst(2, 2, 3, 2, {1, 2}));
  OMEGA_H_CHECK(compare_hst(2, 2, 3, 3, {0, 3}));
  OMEGA_H_CHECK(compare_hst(3, 0, 0, 0, {3, 0}));
  for (Int i = 0; i <= 6; ++i) {
    OMEGA_H_CHECK(compare_hst(3, 1, i, 0, {2, i}));
    OMEGA_H_CHECK(compare_hst(3, 1, i, 1, {3, 0}));
  }
  OMEGA_H_CHECK(compare_hst(3, 2, 0, 0, {1, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 0, 1, {2, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 0, 2, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 0, 3, {2, 1}));
  OMEGA_H_CHECK(compare_hst(3, 2, 1, 0, {2, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 1, 1, {1, 1}));
  OMEGA_H_CHECK(compare_hst(3, 2, 1, 2, {2, 2}));
  OMEGA_H_CHECK(compare_hst(3, 2, 1, 3, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 2, 0, {2, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 2, 1, {1, 2}));
  OMEGA_H_CHECK(compare_hst(3, 2, 2, 2, {2, 3}));
  OMEGA_H_CHECK(compare_hst(3, 2, 2, 3, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 3, 0, {1, 3}));
  OMEGA_H_CHECK(compare_hst(3, 2, 3, 1, {2, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 3, 2, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 3, 3, {2, 4}));
  OMEGA_H_CHECK(compare_hst(3, 2, 4, 0, {1, 4}));
  OMEGA_H_CHECK(compare_hst(3, 2, 4, 1, {2, 1}));
  OMEGA_H_CHECK(compare_hst(3, 2, 4, 2, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 4, 3, {2, 4}));
  OMEGA_H_CHECK(compare_hst(3, 2, 5, 0, {2, 1}));
  OMEGA_H_CHECK(compare_hst(3, 2, 5, 1, {1, 5}));
  OMEGA_H_CHECK(compare_hst(3, 2, 5, 2, {2, 2}));
  OMEGA_H_CHECK(compare_hst(3, 2, 5, 3, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 6, 0, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 6, 1, {2, 2}));
  OMEGA_H_CHECK(compare_hst(3, 2, 6, 2, {1, 6}));
  OMEGA_H_CHECK(compare_hst(3, 2, 6, 3, {2, 3}));
  OMEGA_H_CHECK(compare_hst(3, 2, 7, 0, {2, 4}));
  OMEGA_H_CHECK(compare_hst(3, 2, 7, 1, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 7, 2, {2, 3}));
  OMEGA_H_CHECK(compare_hst(3, 2, 7, 3, {1, 7}));
  OMEGA_H_CHECK(compare_hst(3, 2, 8, 0, {2, 1}));
  OMEGA_H_CHECK(compare_hst(3, 2, 8, 1, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 8, 2, {2, 5}));
  OMEGA_H_CHECK(compare_hst(3, 2, 8, 3, {1, 8}));
  OMEGA_H_CHECK(compare_hst(3, 2, 9, 0, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 9, 1, {2, 2}));
  OMEGA_H_CHECK(compare_hst(3, 2, 9, 2, {1, 9}));
  OMEGA_H_CHECK(compare_hst(3, 2, 9, 3, {2, 5}));
  OMEGA_H_CHECK(compare_hst(3, 2, 10, 0, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 10, 1, {2, 3}));
  OMEGA_H_CHECK(compare_hst(3, 2, 10, 2, {1, 10}));
  OMEGA_H_CHECK(compare_hst(3, 2, 10, 3, {2, 5}));
  OMEGA_H_CHECK(compare_hst(3, 2, 11, 0, {2, 4}));
  OMEGA_H_CHECK(compare_hst(3, 2, 11, 1, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 2, 11, 2, {2, 5}));
  OMEGA_H_CHECK(compare_hst(3, 2, 11, 3, {1, 11}));
  OMEGA_H_CHECK(compare_hst(3, 3, 0, 0, {0, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 0, 1, {1, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 0, 2, {2, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 0, 3, {1, 3}));
  OMEGA_H_CHECK(compare_hst(3, 3, 0, 4, {1, 4}));
  OMEGA_H_CHECK(compare_hst(3, 3, 0, 5, {2, 1}));
  OMEGA_H_CHECK(compare_hst(3, 3, 0, 6, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 0, 7, {2, 4}));
  OMEGA_H_CHECK(compare_hst(3, 3, 1, 0, {1, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 1, 1, {0, 1}));
  OMEGA_H_CHECK(compare_hst(3, 3, 1, 2, {1, 1}));
  OMEGA_H_CHECK(compare_hst(3, 3, 1, 3, {2, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 1, 4, {2, 1}));
  OMEGA_H_CHECK(compare_hst(3, 3, 1, 5, {1, 5}));
  OMEGA_H_CHECK(compare_hst(3, 3, 1, 6, {2, 2}));
  OMEGA_H_CHECK(compare_hst(3, 3, 1, 7, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 2, 0, {2, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 2, 1, {1, 1}));
  OMEGA_H_CHECK(compare_hst(3, 3, 2, 2, {0, 2}));
  OMEGA_H_CHECK(compare_hst(3, 3, 2, 3, {1, 2}));
  OMEGA_H_CHECK(compare_hst(3, 3, 2, 4, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 2, 5, {2, 2}));
  OMEGA_H_CHECK(compare_hst(3, 3, 2, 6, {1, 6}));
  OMEGA_H_CHECK(compare_hst(3, 3, 2, 7, {2, 3}));
  OMEGA_H_CHECK(compare_hst(3, 3, 3, 0, {1, 3}));
  OMEGA_H_CHECK(compare_hst(3, 3, 3, 1, {2, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 3, 2, {1, 2}));
  OMEGA_H_CHECK(compare_hst(3, 3, 3, 3, {0, 3}));
  OMEGA_H_CHECK(compare_hst(3, 3, 3, 4, {2, 4}));
  OMEGA_H_CHECK(compare_hst(3, 3, 3, 5, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 3, 6, {2, 3}));
  OMEGA_H_CHECK(compare_hst(3, 3, 3, 7, {1, 7}));
  OMEGA_H_CHECK(compare_hst(3, 3, 4, 0, {1, 4}));
  OMEGA_H_CHECK(compare_hst(3, 3, 4, 1, {2, 1}));
  OMEGA_H_CHECK(compare_hst(3, 3, 4, 2, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 4, 3, {2, 4}));
  OMEGA_H_CHECK(compare_hst(3, 3, 4, 4, {0, 4}));
  OMEGA_H_CHECK(compare_hst(3, 3, 4, 5, {1, 8}));
  OMEGA_H_CHECK(compare_hst(3, 3, 4, 6, {2, 5}));
  OMEGA_H_CHECK(compare_hst(3, 3, 4, 7, {1, 11}));
  OMEGA_H_CHECK(compare_hst(3, 3, 5, 0, {2, 1}));
  OMEGA_H_CHECK(compare_hst(3, 3, 5, 1, {1, 5}));
  OMEGA_H_CHECK(compare_hst(3, 3, 5, 2, {2, 2}));
  OMEGA_H_CHECK(compare_hst(3, 3, 5, 3, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 5, 4, {1, 8}));
  OMEGA_H_CHECK(compare_hst(3, 3, 5, 5, {0, 5}));
  OMEGA_H_CHECK(compare_hst(3, 3, 5, 6, {1, 9}));
  OMEGA_H_CHECK(compare_hst(3, 3, 5, 7, {2, 5}));
  OMEGA_H_CHECK(compare_hst(3, 3, 6, 0, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 6, 1, {2, 2}));
  OMEGA_H_CHECK(compare_hst(3, 3, 6, 2, {1, 6}));
  OMEGA_H_CHECK(compare_hst(3, 3, 6, 3, {2, 3}));
  OMEGA_H_CHECK(compare_hst(3, 3, 6, 4, {2, 5}));
  OMEGA_H_CHECK(compare_hst(3, 3, 6, 5, {1, 9}));
  OMEGA_H_CHECK(compare_hst(3, 3, 6, 6, {0, 6}));
  OMEGA_H_CHECK(compare_hst(3, 3, 6, 7, {1, 10}));
  OMEGA_H_CHECK(compare_hst(3, 3, 7, 0, {2, 4}));
  OMEGA_H_CHECK(compare_hst(3, 3, 7, 1, {3, 0}));
  OMEGA_H_CHECK(compare_hst(3, 3, 7, 2, {2, 3}));
  OMEGA_H_CHECK(compare_hst(3, 3, 7, 3, {1, 7}));
  OMEGA_H_CHECK(compare_hst(3, 3, 7, 4, {1, 11}));
  OMEGA_H_CHECK(compare_hst(3, 3, 7, 5, {2, 5}));
  OMEGA_H_CHECK(compare_hst(3, 3, 7, 6, {1, 10}));
  OMEGA_H_CHECK(compare_hst(3, 3, 7, 7, {0, 7}));
}

void test_copy_constructor(Library* lib)
{
  auto world = lib->world();
  fprintf(stderr, "before build_box\n");
  auto mesh_a = build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 1, 1, 1,false);
  fprintf(stderr, "after build_box mesh_a at %p\n", &mesh_a);
  auto mesh_b = mesh_a;
  fprintf(stderr, "mesh_b = mesh_a, mesh_b at %p\n", &mesh_b);
  OMEGA_H_CHECK(mesh_a.coords() == mesh_b.coords());
  Write<Real> two_coords_w(mesh_a.coords().size());
  Omega_h::parallel_for(two_coords_w.size(), OMEGA_H_LAMBDA(LO i){ two_coords_w[i]=2.0; });
  Read<Real> two_coords(two_coords_w);
  mesh_b.set_coords(two_coords);
  OMEGA_H_CHECK(!(mesh_a.coords() == mesh_b.coords()));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  test_down_template();
  test_tri_align();
  test_form_uses();
  test_reflect_down();
  test_find_unique();
  test_hilbert();
  test_bbox();
  test_build(&lib);
  test_star(&lib);
  test_dual(&lib);
  test_quality();
  test_inertial_bisect(&lib);
  test_average_field(&lib);
  test_refine_qualities(&lib);
  test_mark_up_down(&lib);
  test_compare_meshes(&lib);
  test_swap2d_topology(&lib);
  test_swap3d_loop(&lib);
  test_element_implied_metric();
  test_recover_hessians(&lib);
  test_sf_scale(&lib);
  test_proximity(&lib);
  test_1d_box(&lib);
  test_hypercube_split_template();
  test_copy_constructor(&lib);
}
