#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_math.hpp"
#include "Omega_h_motion.hpp"
#include "Omega_h_proximity.hpp"
#include "adjacency.hpp"
#include "align.hpp"
#include "bbox.hpp"
#include "derive.hpp"
#include "eigen.hpp"
#include "file.hpp"
#include "hilbert.hpp"
#include "inertia.hpp"
#include "int128.hpp"
#include "internal.hpp"
#include "lie.hpp"
#include "linpart.hpp"
#include "loop.hpp"
#include "mark.hpp"
#include "quality.hpp"
#include "refine_qualities.hpp"
#include "scan.hpp"
#include "size.hpp"
#include "sort.hpp"
#include "swap2d.hpp"
#include "swap3d_choice.hpp"
#include "swap3d_loop.hpp"
#include "transfer_conserve.hpp"
#include "vtk.hpp"
#include "xml.hpp"

#include <sstream>

using namespace Omega_h;

template <Int m, Int n>
static void test_qr_decomp(Matrix<m, n> a) {
  auto qr = factorize_qr_householder(m, n, a);
  auto r = qr.r;
  auto q = identity_matrix<m, n>();
  for (Int j = 0; j < n; ++j) implicit_q_x(m, n, q[j], qr.v);
  CHECK(are_close(a, q * r));
  CHECK(are_close(transpose(q) * q, identity_matrix<n, n>()));
}

static void test_qr_decomps() {
  test_qr_decomp(identity_matrix<3, 3>());
  test_qr_decomp(Matrix<3, 3>({EPSILON, 0, 0, 0, EPSILON, 0, 0, 0, EPSILON}));
  test_qr_decomp(Matrix<3, 3>({12, -51, 4, 6, 167, -68, -4, 24, -41}));
}

static void test_form_ortho_basis() {
  auto n = normalize(vector_3(1, 1, 1));
  auto f = form_ortho_basis(n);
  CHECK(are_close(f[0], n));
  CHECK(are_close(transpose(f) * f, identity_matrix<3, 3>()));
}

static void test_least_squares() {
  Matrix<4, 2> m({1, 1, 1, 2, 1, 3, 1, 4});
  Vector<4> b({6, 5, 7, 10});
  auto x = solve_using_qr(m, b);
  CHECK(are_close(x, vector_2(3.5, 1.4)));
}

static void test_int128() {
  Int128 a(INT64_MAX);
  auto b = a + a;
  b = b + b;
  b = b + b;
  b = b >> 3;
  CHECK(b == a);
}

static void test_repro_sum() {
  Reals a({std::exp2(int(20)), std::exp2(int(-20))});
  Real sum = repro_sum(a);
  CHECK(sum == std::exp2(20) + std::exp2(int(-20)));
}

static void test_cubic(Few<Real, 3> coeffs, Int nroots_wanted,
    Few<Real, 3> roots_wanted, Few<Int, 3> mults_wanted) {
  auto roots = find_polynomial_roots(coeffs);
  CHECK(roots.n == nroots_wanted);
  for (Int i = 0; i < roots.n; ++i) {
    CHECK(roots.mults[i] == mults_wanted[i]);
    CHECK(are_close(roots.values[i], roots_wanted[i]));
  }
}

static void test_cubic() {
  test_cubic({0, 0, 0}, 1, Few<Real, 3>({0}), Few<Int, 3>({3}));
  test_cubic({1., -3. / 2., -3. / 2.}, 3, Few<Real, 3>({2, -1, .5}),
      Few<Int, 3>({1, 1, 1}));
  test_cubic({2., -3., 0}, 2, Few<Real, 3>({-2, 1}), Few<Int, 3>({1, 2}));
  test_cubic({-8, -6, 3}, 3, Few<Real, 3>({2, -4, -1}), Few<Int, 3>({1, 1, 1}));
}

template <Int dim>
static void test_eigen(
    Matrix<dim, dim> m, Matrix<dim, dim> q_expect, Vector<dim> l_expect) {
  auto ed = decompose_eigen(m);
  auto q = ed.q;
  auto l = ed.l;
  CHECK(are_close(q, q_expect));
  CHECK(are_close(l, l_expect));
}

template <Int dim>
static void test_eigen(Matrix<dim, dim> m, Vector<dim> l_expect) {
  auto ed = decompose_eigen(m);
  auto q = ed.q;
  auto l = ed.l;
  CHECK(are_close(l, l_expect, 1e-8, 1e-8));
  CHECK(are_close(m, compose_eigen(q, l)));
}

static void test_eigen_cubic_ortho(Matrix<3, 3> m, Vector<3> l_expect) {
  auto ed = decompose_eigen(m);
  auto q = ed.q;
  auto l = ed.l;
  CHECK(are_close(transpose(q) * q, identity_matrix<3, 3>(), 1e-8, 1e-8));
  CHECK(are_close(l, l_expect, 1e-8, 1e-8));
  CHECK(are_close(m, compose_ortho(q, l), 1e-8, 1e-8));
}

static void test_eigen_metric(Vector<3> h) {
  auto q =
      rotate(PI / 4., vector_3(0, 0, 1)) * rotate(PI / 4., vector_3(0, 1, 0));
  CHECK(are_close(transpose(q) * q, identity_matrix<3, 3>()));
  auto l = metric_eigenvalues(h);
  auto a = compose_ortho(q, l);
  test_eigen_cubic_ortho(a, l);
}

static void test_eigen_quadratic() {
  test_eigen(identity_matrix<2, 2>(), identity_matrix<2, 2>(), vector_2(1, 1));
  test_eigen(zero_matrix<2, 2>(), identity_matrix<2, 2>(), vector_2(0, 0));
  test_eigen(matrix_2x2(8.67958, -14.0234, -1.04985, 2.25873),
      matrix_2x2(9.9192948778227130e-01, 8.6289280817702185e-01,
          -1.2679073810022995e-01, 5.0538698202107812e-01),
      vector_2(1.0472083659357935e+01, 4.6622634064206342e-01));
}

static void test_eigen_cubic() {
  test_eigen(
      identity_matrix<3, 3>(), identity_matrix<3, 3>(), vector_3(1, 1, 1));
  test_eigen(zero_matrix<3, 3>(), identity_matrix<3, 3>(), vector_3(0, 0, 0));
  test_eigen(matrix_3x3(-1, 3, -1, -3, 5, -1, -3, 3, 1), vector_3(1, 2, 2));
  /* the lengths have to be ordered so that
     if two of them are the same they should
     appear at the end */
  test_eigen_metric(vector_3(1e+3, 1, 1));
  test_eigen_metric(vector_3(1, 1e+3, 1e+3));
  test_eigen_metric(vector_3(1e-3, 1, 1));
  test_eigen_metric(vector_3(1, 1e-3, 1e-3));
  test_eigen_metric(vector_3(1e-6, 1e-3, 1e-3));
}

static void test_intersect_ortho_metrics(
    Vector<3> h1, Vector<3> h2, Vector<3> hi_expect) {
  auto q =
      rotate(PI / 4., vector_3(0, 0, 1)) * rotate(PI / 4., vector_3(0, 1, 0));
  auto m1 = compose_metric(q, h1);
  auto m2 = compose_metric(q, h2);
  auto mi = intersect_metrics(m1, m2);
  /* if we decompose it, the eigenvectors may
     get re-ordered. */
  for (Int i = 0; i < 3; ++i) {
    CHECK(are_close(metric_desired_length(mi, q[i]), hi_expect[i], 1e-3));
  }
}

static void test_intersect_subset_metrics() {
  auto h1 = vector_2(1, 2);
  auto r1 = identity_matrix<2, 2>();
  auto h2 = vector_2(2, 3);
  auto r2 = rotate(PI / 4);
  auto m1 = compose_metric(r1, h1);
  auto m2 = compose_metric(r2, h2);
  CHECK(are_close(intersect_metrics(m2, m1), m1));
  CHECK(are_close(intersect_metrics(m1, m2), m1));
}

static void test_intersect_metrics() {
  test_intersect_ortho_metrics(
      vector_3(0.5, 1, 1), vector_3(1, 0.5, 1), vector_3(0.5, 0.5, 1));
  test_intersect_ortho_metrics(
      vector_3(1e-2, 1, 1), vector_3(1, 1, 1e-2), vector_3(1e-2, 1, 1e-2));
  test_intersect_ortho_metrics(vector_3(1e-2, 1e-2, 1), vector_3(1, 1, 1e-2),
      vector_3(1e-2, 1e-2, 1e-2));
  test_intersect_ortho_metrics(vector_3(1e-5, 1e-3, 1e-3),
      vector_3(1e-3, 1e-3, 1e-5), vector_3(1e-5, 1e-3, 1e-5));
  test_intersect_subset_metrics();
}

static void test_sort() {
  {
    LOs a({0, 1});
    LOs perm = sort_by_keys(a);
    CHECK(perm == LOs({0, 1}));
  }
  {
    LOs a({0, 2, 0, 1});
    LOs perm = sort_by_keys(a, 2);
    CHECK(perm == LOs({1, 0}));
  }
  {
    LOs a({0, 2, 1, 1});
    LOs perm = sort_by_keys(a, 2);
    CHECK(perm == LOs({0, 1}));
  }
  {
    LOs a({1, 2, 3, 1, 2, 2, 3, 0, 0});
    LOs perm = sort_by_keys(a, 3);
    CHECK(perm == LOs({1, 0, 2}));
  }
}

static void test_scan() {
  {
    LOs scanned = offset_scan(LOs(3, 1));
    CHECK(scanned == Read<LO>(4, 0, 1));
  }
  {
    LOs scanned = offset_scan(Read<I8>(3, 1));
    CHECK(scanned == Read<LO>(4, 0, 1));
  }
}

static void test_fan_and_funnel() {
  CHECK(invert_funnel(LOs({0, 0, 1, 1, 2, 2}), 3) == LOs({0, 2, 4, 6}));
  CHECK(invert_fan(LOs({0, 2, 4, 6})) == LOs({0, 0, 1, 1, 2, 2}));
  CHECK(invert_funnel(LOs({0, 0, 0, 2, 2, 2}), 3) == LOs({0, 3, 3, 6}));
  CHECK(invert_fan(LOs({0, 3, 3, 6})) == LOs({0, 0, 0, 2, 2, 2}));
  CHECK(invert_funnel(LOs({0, 0, 0, 0, 0, 0}), 3) == LOs({0, 6, 6, 6}));
  CHECK(invert_fan(LOs({0, 6, 6, 6})) == LOs({0, 0, 0, 0, 0, 0}));
  CHECK(invert_funnel(LOs({2, 2, 2, 2, 2, 2}), 3) == LOs({0, 0, 0, 6}));
  CHECK(invert_fan(LOs({0, 0, 0, 6})) == LOs({2, 2, 2, 2, 2, 2}));
}

static void test_permute() {
  Reals data({0.1, 0.2, 0.3, 0.4});
  LOs perm({3, 2, 1, 0});
  Reals permuted = unmap(perm, data, 1);
  CHECK(permuted == Reals({0.4, 0.3, 0.2, 0.1}));
  Reals back = permute(permuted, perm, 1);
  CHECK(back == data);
}

// these tests can have degree at most 1
// because map::invert doesn't have to be
// deterministic in local ordering
static void test_invert_map(Graph (*invert_funcptr)(LOs a2b, LO nb)) {
  {
    LOs hl2l({});
    auto l2hl = invert_funcptr(hl2l, 4);
    CHECK(l2hl.a2ab == LOs(5, 0));
    CHECK(l2hl.ab2b == LOs({}));
  }
  {
    LOs hl2l({0, 1, 2, 3});
    auto l2hl = invert_funcptr(hl2l, 4);
    CHECK(l2hl.a2ab == LOs(5, 0, 1));
    CHECK(l2hl.ab2b == LOs(4, 0, 1));
  }
}

static void test_invert_map() {
  test_invert_map(invert_map_by_sorting);
  test_invert_map(invert_map_by_atomics);
}

static void test_invert_adj() {
  Adj tris2verts(LOs({0, 1, 2, 2, 3, 0}));
  Read<GO> tri_globals({0, 1});
  Adj verts2tris = invert_adj(tris2verts, 3, 4, tri_globals);
  CHECK(verts2tris.a2ab == offset_scan(LOs({2, 1, 2, 1})));
  CHECK(verts2tris.ab2b == LOs({0, 1, 0, 0, 1, 1}));
  CHECK(verts2tris.codes ==
        Read<I8>({make_code(0, 0, 0), make_code(0, 0, 2), make_code(0, 0, 1),
            make_code(0, 0, 2), make_code(0, 0, 0), make_code(0, 0, 1)}));
}

static bool same_adj(Int a[], Int b[]) {
  for (Int i = 0; i < 3; ++i)
    if (a[i] != b[i]) return false;
  return true;
}

static void test_tri_align() {
  Int ident[3] = {0, 1, 2};
  Int out[3];
  Int out2[3];
  /* check that flipping and rotating do what we want */
  {
    align_adj<3, Int>(make_code(true, 0, 0), ident, out);
    Int expect[3] = {0, 2, 1};
    CHECK(same_adj(out, expect));
  }
  {
    align_adj<3>(make_code(false, 1, 0), ident, out);
    Int expect[3] = {2, 0, 1};
    CHECK(same_adj(out, expect));
  }
  {
    align_adj<3>(make_code(false, 2, 0), ident, out);
    Int expect[3] = {1, 2, 0};
    CHECK(same_adj(out, expect));
  }
  /* check that compound_alignments does its job */
  for (I8 rot1 = 0; rot1 < 3; ++rot1)
    for (I8 flip1 = 0; flip1 < 2; ++flip1)
      for (I8 rot2 = 0; rot2 < 3; ++rot2)
        for (I8 flip2 = 0; flip2 < 2; ++flip2) {
          I8 code1 = make_code(flip1, rot1, 0);
          I8 code2 = make_code(flip2, rot2, 0);
          align_adj<3>(code1, ident, out);
          align_adj<3>(code2, out, out2);
          Int out3[3];
          I8 code3 = compound_alignments<3>(code1, code2);
          align_adj<3>(code3, ident, out3);
          CHECK(same_adj(out2, out3));
        }
}

static void test_form_uses() {
  CHECK(form_uses(LOs({0, 1, 2}), 2, 1) == LOs({0, 1, 1, 2, 2, 0}));
  CHECK(form_uses(LOs({0, 1, 2, 3}), 3, 1) ==
        LOs({0, 1, 1, 2, 2, 0, 0, 3, 1, 3, 2, 3}));
  CHECK(form_uses(LOs({0, 1, 2, 3}), 3, 2) ==
        LOs({0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3}));
}

static void test_reflect_down() {
  Adj a;
  a = reflect_down(LOs({}), LOs({}), 0, 2, 1);
  CHECK(a.ab2b == LOs({}));
  CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({}), LOs({}), 0, 3, 1);
  CHECK(a.ab2b == LOs({}));
  CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({}), LOs({}), 0, 3, 2);
  CHECK(a.ab2b == LOs({}));
  CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({0, 1, 2}), LOs({0, 1, 1, 2, 2, 0}), 3, 2, 1);
  CHECK(a.ab2b == LOs({0, 1, 2}));
  CHECK(a.codes == Read<I8>({0, 0, 0}));
  a = reflect_down(
      LOs({0, 1, 2, 3}), LOs({0, 1, 1, 2, 2, 0, 0, 3, 1, 3, 2, 3}), 4, 3, 1);
  CHECK(a.ab2b == LOs({0, 1, 2, 3, 4, 5}));
  CHECK(a.codes == Read<I8>({0, 0, 0, 0, 0, 0}));
  a = reflect_down(
      LOs({0, 1, 2, 3}), LOs({0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3}), 4, 3, 2);
  CHECK(a.ab2b == LOs({0, 1, 2, 3}));
  CHECK(a.codes == Read<I8>({0, 0, 0, 0}));
  a = reflect_down(
      LOs({0, 1, 2, 3}), LOs({0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0}), 4, 3, 2);
  CHECK(a.ab2b == LOs({0, 1, 2, 3}));
  CHECK(a.codes == Read<I8>(4, make_code(true, 0, 0)));
  a = reflect_down(
      LOs({0, 1, 2, 2, 3, 0}), LOs({0, 1, 1, 2, 2, 3, 3, 0, 0, 2}), 4, 2, 1);
  CHECK(a.ab2b == LOs({0, 1, 4, 2, 3, 4}));
}

static void test_find_unique() {
  CHECK(find_unique(LOs({}), 2, 1) == LOs({}));
  CHECK(find_unique(LOs({}), 3, 1) == LOs({}));
  CHECK(find_unique(LOs({}), 3, 2) == LOs({}));
  CHECK(find_unique(LOs({0, 1, 2, 2, 3, 0}), 2, 1) ==
        LOs({0, 1, 0, 2, 3, 0, 1, 2, 2, 3}));
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
  CHECK(stream.str() == expected);
  hilbert::coord_t Y[3];
  hilbert::untranspose(X, Y, 5, 3);
  std::stringstream stream2;
  stream2 << "Hilbert integer = " << (Y[0] >> 4 & 1) << (Y[0] >> 3 & 1)
          << (Y[0] >> 2 & 1) << (Y[0] >> 1 & 1) << (Y[0] >> 0 & 1)
          << (Y[1] >> 4 & 1) << (Y[1] >> 3 & 1) << (Y[1] >> 2 & 1)
          << (Y[1] >> 1 & 1) << (Y[1] >> 0 & 1) << (Y[2] >> 4 & 1)
          << (Y[2] >> 3 & 1) << (Y[2] >> 2 & 1) << (Y[2] >> 1 & 1)
          << (Y[2] >> 0 & 1) << " = 7865 check";
  CHECK(stream2.str() == expected);
}

static void test_bbox() {
  CHECK(are_close(BBox<2>(vector_2(-3, -3), vector_2(3, 3)),
      find_bounding_box<2>(Reals({0, -3, 3, 0, 0, 3, -3, 0}))));
  CHECK(are_close(BBox<3>(vector_3(-3, -3, -3), vector_3(3, 3, 3)),
      find_bounding_box<3>(
          Reals({0, -3, 0, 3, 0, 0, 0, 3, 0, -3, 0, 0, 0, 0, -3, 0, 0, 3}))));
}

static void test_build_from_elems2verts(Library* lib) {
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 2, LOs({0, 1, 2}), 3);
    CHECK(mesh.ask_down(2, 0).ab2b == LOs({0, 1, 2}));
    CHECK(mesh.ask_down(2, 1).ab2b == LOs({0, 2, 1}));
    CHECK(mesh.ask_down(1, 0).ab2b == LOs({0, 1, 2, 0, 1, 2}));
  }
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 3, LOs({0, 1, 2, 3}), 4);
    CHECK(mesh.ask_down(3, 0).ab2b == LOs({0, 1, 2, 3}));
  }
}

static void test_star(Library* lib) {
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 2, LOs({0, 1, 2}), 3);
    Adj v2v = mesh.ask_star(VERT);
    CHECK(v2v.a2ab == LOs(4, 0, 2));
    CHECK(v2v.ab2b == LOs({1, 2, 0, 2, 0, 1}));
    Adj e2e = mesh.ask_star(EDGE);
    CHECK(e2e.a2ab == LOs(4, 0, 2));
    CHECK(e2e.ab2b == LOs({2, 1, 0, 2, 1, 0}));
  }
  {
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 3, LOs({0, 1, 2, 3}), 4);
    Adj v2v = mesh.ask_star(VERT);
    CHECK(v2v.a2ab == LOs(5, 0, 3));
    CHECK(v2v.ab2b == LOs({1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2}));
    Adj e2e = mesh.ask_star(EDGE);
    CHECK(e2e.a2ab == LOs(7, 0, 5));
    CHECK(e2e.ab2b == LOs({1, 3, 4, 2, 5, 3, 0, 2, 5, 4, 0, 4, 5, 1, 3, 0, 1, 5,
                          4, 2, 2, 0, 3, 5, 1, 1, 2, 4, 3, 0}));
  }
}

static void test_injective_map() {
  LOs primes2ints({2, 3, 5, 7});
  LOs ints2primes = invert_injective_map(primes2ints, 8);
  CHECK(ints2primes == LOs({-1, -1, 0, 1, -1, 2, -1, 3}));
}

static void test_dual(Library* lib) {
  Mesh mesh(lib);
  build_from_elems2verts(&mesh, 2, LOs({0, 1, 2, 2, 3, 0}), 4);
  auto t2t = mesh.ask_dual();
  auto t2tt = t2t.a2ab;
  auto tt2t = t2t.ab2b;
  CHECK(t2tt == offset_scan(LOs({1, 1})));
  CHECK(tt2t == LOs({1, 0}));
}

static void test_quality() {
  Few<Vector<2>, 3> perfect_tri(
      {vector_2(1, 0), vector_2(0, sqrt(3.0)), vector_2(-1, 0)});
  Few<Vector<3>, 4> perfect_tet(
      {vector_3(1, 0, -1.0 / sqrt(2.0)), vector_3(-1, 0, -1.0 / sqrt(2.0)),
          vector_3(0, -1, 1.0 / sqrt(2.0)), vector_3(0, 1, 1.0 / sqrt(2.0))});
  Few<Vector<2>, 3> flat_tri({vector_2(1, 0), vector_2(0, 0), vector_2(-1, 0)});
  Few<Vector<3>, 4> flat_tet({vector_3(1, 0, 0), vector_3(-1, 0, 0),
      vector_3(0, -1, 0), vector_3(0, 1, 0)});
  Few<Vector<2>, 3> inv_tri(
      {vector_2(1, 0), vector_2(-1, 0), vector_2(0, sqrt(3.0))});
  Few<Vector<3>, 4> inv_tet(
      {vector_3(1, 0, -1.0 / sqrt(2.0)), vector_3(-1, 0, -1.0 / sqrt(2.0)),
          vector_3(0, 1, 1.0 / sqrt(2.0)), vector_3(0, -1, 1.0 / sqrt(2.0))});
  Matrix<2, 2> id_metric_2 = identity_matrix<2, 2>();
  Matrix<3, 3> id_metric_3 = identity_matrix<3, 3>();
  Matrix<2, 2> x_metric_2 =
      compose_metric(identity_matrix<2, 2>(), vector_2(1, 0.5));
  Matrix<3, 3> x_metric_3 =
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
  CHECK(are_close(real_element_quality(perfect_tri), 1.0));
  CHECK(are_close(real_element_quality(perfect_tet), 1.0));
  CHECK(are_close(real_element_quality(flat_tri), 0.0));
  CHECK(are_close(real_element_quality(flat_tet), 0.0));
  CHECK(real_element_quality(inv_tri) < 0.0);
  CHECK(real_element_quality(inv_tet) < 0.0);
  CHECK(are_close(metric_element_quality(perfect_tri, id_metric_2), 1.0));
  CHECK(are_close(metric_element_quality(perfect_tet, id_metric_3), 1.0));
  CHECK(are_close(metric_element_quality(flat_tri, id_metric_2), 0.0));
  CHECK(are_close(metric_element_quality(flat_tet, id_metric_3), 0.0));
  CHECK(metric_element_quality(inv_tri, id_metric_2) < 0.0);
  CHECK(metric_element_quality(inv_tet, id_metric_3) < 0.0);
  CHECK(are_close(metric_element_quality(x_tri, x_metric_2), 1.0));
  CHECK(are_close(metric_element_quality(x_tet, x_metric_3), 1.0));
}

static void test_file_components() {
  using namespace binary;
  std::stringstream stream;
  std::string s = "foo";
  LO n = 10;
#ifdef OMEGA_H_USE_ZLIB
  bool is_compressed = true;
#else
  bool is_compressed = false;
#endif
  I8 a = 2;
  write_value(stream, a);
  I32 b = 42 * 1000 * 1000;
  write_value(stream, b);
  I64 c = I64(42) * 1000 * 1000 * 1000;
  write_value(stream, c);
  Real d = 4.2;
  write_value(stream, d);
  Read<I8> aa(n, 0, a);
  write_array(stream, aa);
  Read<I32> ab(n, 0, b);
  write_array(stream, ab);
  Read<I64> ac(n, 0, c);
  write_array(stream, ac);
  Read<Real> ad(n, 0, d);
  write_array(stream, ad);
  write(stream, s);
  I8 a2;
  read_value(stream, a2);
  CHECK(a == a2);
  I32 b2;
  read_value(stream, b2);
  CHECK(b == b2);
  I64 c2;
  read_value(stream, c2);
  CHECK(c == c2);
  Real d2;
  read_value(stream, d2);
  CHECK(d == d2);
  Read<I8> aa2;
  read_array(stream, aa2, is_compressed);
  CHECK(aa2 == aa);
  Read<I32> ab2;
  read_array(stream, ab2, is_compressed);
  CHECK(ab2 == ab);
  Read<I64> ac2;
  read_array(stream, ac2, is_compressed);
  CHECK(ac2 == ac);
  Read<Real> ad2;
  read_array(stream, ad2, is_compressed);
  CHECK(ad2 == ad);
  std::string s2;
  read(stream, s2);
  CHECK(s == s2);
}

static void test_linpart() {
  GO total = 7;
  I32 comm_size = 2;
  CHECK(linear_partition_size(total, comm_size, 0) == 4);
  CHECK(linear_partition_size(total, comm_size, 1) == 3);
  Read<GO> globals({6, 5, 4, 3, 2, 1, 0});
  auto remotes = globals_to_linear_owners(globals, total, comm_size);
  CHECK(remotes.ranks == Read<I32>({1, 1, 1, 0, 0, 0, 0}));
  CHECK(remotes.idxs == Read<I32>({2, 1, 0, 3, 2, 1, 0}));
}

static void test_expand() {
  auto fan = offset_scan(LOs({2, 1, 3}));
  Reals data({2.2, 3.14, 42.0});
  CHECK(expand(data, fan, 1) == Reals({2.2, 2.2, 3.14, 42.0, 42.0, 42.0}));
}

static void test_inertial_bisect(Library* lib) {
  Reals coords({2, 1, 0, 2, -1, 0, -2, 1, 0, -2, -1, 0});
  Reals masses(4, 1);
  auto self = lib->self();
  Real tolerance = 0.0;
  Vector<3> axis;
  auto marked = inertia::mark_bisection(self, coords, masses, tolerance, axis);
  CHECK(marked == Read<I8>({1, 1, 0, 0}));
  marked = inertia::mark_bisection_given_axis(
      self, coords, masses, tolerance, vector_3(0, 1, 0));
  CHECK(marked == Read<I8>({1, 0, 1, 0}));
}

static void test_average_field(Library* lib) {
  Mesh mesh(lib);
  build_box(&mesh, 1, 1, 0, 1, 1, 0);
  Reals v2x({2, 1, 3, 2});
  auto e2x = average_field(&mesh, 2, LOs({0, 1}), 1, v2x);
  CHECK(are_close(e2x, Reals({5. / 3., 7. / 3.})));
}

template <Int n>
static void test_positivize(Vector<n> pos) {
  auto neg = pos * -1.0;
  CHECK(are_close(positivize(pos), pos));
  CHECK(are_close(positivize(neg), pos));
}

static void test_positivize() {
  test_positivize(vector_3(1, 1, 1));
  test_positivize(vector_3(1, -1, 1));
  test_positivize(vector_2(-1, 1));
  test_positivize(vector_2(1, 1));
}

static void test_edge_length() {
  CHECK(are_close(1., edge_length(1., 1.)));
  CHECK(edge_length(1., 2.) > 1.);
  CHECK(edge_length(1., 2.) < 1.5);
}

static void test_refine_qualities(Library* lib) {
  Mesh mesh(lib);
  build_box(&mesh, 1, 1, 0, 1, 1, 0);
  LOs candidates(mesh.nedges(), 0, 1);
  auto quals = refine_qualities(&mesh, candidates);
  CHECK(are_close(
      quals, Reals({0.494872, 0.494872, 0.866025, 0.494872, 0.494872}), 1e-4));
  mesh.add_tag(VERT, "metric", symm_dofs(2), OMEGA_H_METRIC, OMEGA_H_DO_OUTPUT,
      repeat_symm(mesh.nverts(), identity_matrix<2, 2>()));
  auto quals2 = refine_qualities(&mesh, candidates);
  CHECK(are_close(quals2, quals));
}

static void test_mark_up_down(Library* lib) {
  Mesh mesh(lib);
  build_box(&mesh, 1, 1, 0, 1, 1, 0);
  CHECK(
      mark_down(&mesh, TRI, VERT, Read<I8>({1, 0})) == Read<I8>({1, 1, 0, 1}));
  CHECK(mark_up(&mesh, VERT, TRI, Read<I8>({0, 1, 0, 0})) == Read<I8>({1, 0}));
}

static void test_compare_meshes(Library* lib) {
  Mesh a(lib);
  build_box(&a, 1, 1, 0, 4, 4, 0);
  CHECK(a == a);
  Mesh b = a;
  b.reorder();
  CHECK(a == b);
  b.add_tag<I8>(VERT, "foo", 1, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT,
      Read<I8>(b.nverts(), 1));
  CHECK(!(a == b));
}

static void test_swap2d_topology(Library* lib) {
  Mesh mesh(lib);
  build_box(&mesh, 1, 1, 0, 1, 1, 0);
  HostFew<LOs, 3> keys2prods;
  HostFew<LOs, 3> prod_verts2verts;
  auto keys2edges = LOs({2});
  swap2d_topology(&mesh, keys2edges, &keys2prods, &prod_verts2verts);
  CHECK(prod_verts2verts[EDGE] == LOs({2, 1}));
  CHECK(prod_verts2verts[TRI] == LOs({3, 2, 1, 0, 1, 2}));
  CHECK(keys2prods[EDGE] == offset_scan(LOs({1})));
  CHECK(keys2prods[TRI] == offset_scan(LOs({2})));
}

static void test_swap3d_loop(Library* lib) {
  Mesh mesh(lib);
  build_box(&mesh, 1, 1, 1, 1, 1, 1);
  auto edges2tets = mesh.ask_up(EDGE, TET);
  auto edges2edge_tets = edges2tets.a2ab;
  auto edge_tets2tets = edges2tets.ab2b;
  auto edge_tet_codes = edges2tets.codes;
  auto edge_verts2verts = mesh.ask_verts_of(EDGE);
  auto tet_verts2verts = mesh.ask_verts_of(TET);
  auto f = LAMBDA(LO foo) {
    (void)foo;
    LO edge = 6;
    auto loop = swap3d::find_loop(edges2edge_tets, edge_tets2tets,
        edge_tet_codes, edge_verts2verts, tet_verts2verts, edge);
    CHECK(loop.eev2v[0] == 7);
    CHECK(loop.eev2v[1] == 0);
    CHECK(loop.size == 6);
    LO const expect[6] = {2, 3, 1, 5, 4, 6};
    for (Int i = 0; i < loop.size; ++i)
      CHECK(loop.loop_verts2verts[i] == expect[i]);
  };
  parallel_for(LO(1), f);
}

static void build_empty_mesh(Mesh* mesh, Int dim) {
  build_from_elems_and_coords(mesh, dim, LOs({}), Reals({}));
}

static void test_file(Library* lib, Mesh* mesh0) {
  std::stringstream stream;
  binary::write(stream, mesh0);
  Mesh mesh1(lib);
  mesh1.set_comm(lib->self());
  binary::read(stream, &mesh1, binary::latest_version);
  mesh1.set_comm(lib->world());
  auto opts = MeshCompareOpts::init(mesh0, VarCompareOpts::zero_tolerance());
  compare_meshes(mesh0, &mesh1, opts, true, true);
  CHECK(*mesh0 == mesh1);
}

static void test_file(Library* lib) {
  {
    Mesh mesh0(lib);
    build_box(&mesh0, 1, 1, 1, 1, 1, 1);
    test_file(lib, &mesh0);
  }
  {
    Mesh mesh0(lib);
    build_empty_mesh(&mesh0, 3);
    test_file(lib, &mesh0);
  }
}

static void test_xml() {
  xml::Tag tag;
  CHECK(!xml::parse_tag("AQAAAAAAAADABg", &tag));
  CHECK(!xml::parse_tag("   <Foo bar=\"qu", &tag));
  CHECK(!xml::parse_tag("   <Foo bar=", &tag));
  CHECK(xml::parse_tag("   <Foo bar=\"quux\"   >", &tag));
  CHECK(tag.elem_name == "Foo");
  CHECK(tag.attribs["bar"] == "quux");
  CHECK(tag.type == xml::Tag::START);
  CHECK(xml::parse_tag("   <Elem att=\"val\"  answer=\"42\" />", &tag));
  CHECK(tag.elem_name == "Elem");
  CHECK(tag.attribs["att"] == "val");
  CHECK(tag.attribs["answer"] == "42");
  CHECK(tag.type == xml::Tag::SELF_CLOSING);
  CHECK(xml::parse_tag("</Foo>", &tag));
  CHECK(tag.elem_name == "Foo");
  CHECK(tag.type == xml::Tag::END);
}

static void test_read_vtu(Mesh* mesh0) {
  std::stringstream stream;
  vtk::write_vtu(stream, mesh0, mesh0->dim());
  Mesh mesh1(mesh0->library());
  vtk::read_vtu(stream, mesh0->comm(), &mesh1);
  auto opts = MeshCompareOpts::init(mesh0, VarCompareOpts::zero_tolerance());
  CHECK(OMEGA_H_SAME == compare_meshes(mesh0, &mesh1, opts, true, false));
}

static void test_read_vtu(Library* lib) {
  Mesh mesh0(lib);
  build_box(&mesh0, 1, 1, 1, 1, 1, 1);
  test_read_vtu(&mesh0);
}

static void test_interpolate_metrics() {
  auto a = repeat_symm(
      4, compose_metric(identity_matrix<2, 2>(), vector_2(1.0 / 100.0, 1.0)));
  auto b = repeat_symm(
      4, compose_metric(identity_matrix<2, 2>(), vector_2(1.0, 1.0)));
  auto c = interpolate_between_metrics(2, a, b, 0.0);
  CHECK(are_close(a, c));
  c = interpolate_between_metrics(2, a, b, 1.0);
  CHECK(are_close(b, c));
}

static void test_interpolate_isos() {
  auto a = Reals({0.1, 0.2, 0.3});
  auto b = Reals({0.2, 0.4, 0.6});
  auto c = interpolate_between_isos(a, b, 0.0);
  CHECK(are_close(a, c));
  c = interpolate_between_isos(a, b, 1.0);
  CHECK(are_close(b, c));
  c = interpolate_between_isos(a, b, 0.5);
  auto f = LAMBDA(LO i) {
    CHECK(a[i] <= c[i]);
    CHECK(c[i] <= b[i]);
  };
  parallel_for(c.size(), f);
}

static void test_element_implied_metric() {
  /* perfect tri with edge lengths = 2 */
  Few<Vector<2>, 3> perfect_tri(
      {vector_2(1, 0), vector_2(0, sqrt(3.0)), vector_2(-1, 0)});
  auto afm = element_implied_metric(perfect_tri);
  auto bfm = compose_metric(identity_matrix<2, 2>(), vector_2(2, 2));
  CHECK(are_close(afm, bfm));
  /* perfect tet with edge lengths = 2 */
  Few<Vector<3>, 4> perfect_tet(
      {vector_3(1, 0, -1.0 / sqrt(2.0)), vector_3(-1, 0, -1.0 / sqrt(2.0)),
          vector_3(0, -1, 1.0 / sqrt(2.0)), vector_3(0, 1, 1.0 / sqrt(2.0))});
  auto arm = element_implied_metric(perfect_tet);
  auto brm = compose_metric(identity_matrix<3, 3>(), vector_3(2, 2, 2));
  CHECK(are_close(arm, brm));
}

template <Int dim>
static void test_recover_hessians_dim(Library* lib) {
  Mesh mesh(lib);
  Int one_if_3d = ((dim == 3) ? 1 : 0);
  build_box(&mesh, 1, 1, one_if_3d, 4, 4, 4 * one_if_3d);
  classify_by_angles(&mesh, Omega_h::PI / 4);
  auto u_w = Write<Real>(mesh.nverts());
  auto coords = mesh.coords();
  // attach a field = x^2 + y^2 (+ z^2)
  auto f = LAMBDA(LO v) {
    auto x = get_vector<dim>(coords, v);
    u_w[v] = norm_squared(x);
  };
  parallel_for(mesh.nverts(), f);
  auto u = Omega_h::Reals(u_w);
  mesh.add_tag(
      Omega_h::VERT, "u", 1, OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT, u);
  auto hess = recover_hessians(&mesh, u);
  // its second derivative is exactly 2dx + 2dy,
  // and both recovery steps are linear so the current
  // algorithm should get an exact answer
  Vector<dim> dv;
  for (Int i = 0; i < dim; ++i) dv[i] = 2;
  auto expected_hess = repeat_symm(mesh.nverts(), diagonal(dv));
  CHECK(are_close(hess, expected_hess));
}

static void test_recover_hessians(Library* lib) {
  test_recover_hessians_dim<2>(lib);
  test_recover_hessians_dim<3>(lib);
}

template <Int dim>
static void test_sf_scale_dim(Library* lib) {
  Mesh mesh(lib);
  Int one_if_3d = ((dim == 3) ? 1 : 0);
  build_box(&mesh, 1, 1, one_if_3d, 4, 4, 4 * one_if_3d);
  classify_by_angles(&mesh, Omega_h::PI / 4);
  auto target_nelems = mesh.nelems();
  {
    auto size = Omega_h::find_implied_size(&mesh);
    auto size_scal = size_scalar_for_nelems(&mesh, size, target_nelems);
    CHECK(are_close(size_scal, 1.));
  }
  {
    auto metric = Omega_h::find_implied_metric(&mesh);
    auto metric_scal = metric_scalar_for_nelems(&mesh, metric, target_nelems);
    if (dim != 3) CHECK(are_close(metric_scal, 1.));
  }
}

static void test_sf_scale(Library* lib) {
  test_sf_scale_dim<2>(lib);
  test_sf_scale_dim<3>(lib);
}

static void test_categorize_graph() {
  auto g = Graph(LOs({0, 4, 8}), LOs({0, 1, 2, 3, 4, 5, 6, 7}));
  auto b_categories = Read<I32>({8, 8, 42, 8, 42, 42, 42, 42});
  auto g8 = Graph(LOs({0, 3, 3}), LOs({0, 1, 3}));
  auto g42 = Graph(LOs({0, 1, 5}), LOs({2, 4, 5, 6, 7}));
  auto result = categorize_graph(g, b_categories);
  CHECK(result.size() == 2);
  auto tmp = result[8];
  CHECK(result[8] == g8);
  CHECK(result[42] == g42);
}

static void test_circumcenter() {
  Few<Vector<3>, 3> right_tri(
      {vector_3(0, 0, 0), vector_3(1, 0, 0), vector_3(0, 1, 0)});
  auto v0 = get_circumcenter_vector(simplex_basis<3, 2>(right_tri));
  CHECK(are_close(v0, vector_3(0.5, 0.5, 0)));
  Few<Vector<3>, 3> equal_tri(
      {vector_3(0, sqrt(3), 0), vector_3(-1, 0, 0), vector_3(1, 0, 0)});
  auto v1 = get_circumcenter_vector(simplex_basis<3, 2>(equal_tri));
  CHECK(are_close(v1, vector_3(0, -sqrt(3) * 2.0 / 3.0, 0)));
}

static void test_lie() {
  auto a = identity_matrix<3, 3>();
  auto log_a = log_glp(a);
  auto a2 = exp_glp(log_a);
  CHECK(are_close(a2, a));
}

static void test_proximity(Library* lib) {
  {  // triangle with one bridge
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 2, LOs({0, 1, 2}), 3);
    mesh.add_tag(VERT, "coordinates", 2, OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals({0, 0, 1, 0, 0, 1}));
    auto isos = get_pad_isos(&mesh, 2, 1.0, 42.0, Read<I8>({0, 1, 0}));
    CHECK(isos == Reals({42.0}));
  }
  {  // triangle off-center
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 2, LOs({0, 1, 2}), 3);
    mesh.add_tag(VERT, "coordinates", 2, OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals({0, 0, 1, 1, 1, 2}));
    auto isos = get_pad_isos(&mesh, 2, 1.0, 42.0, Read<I8>({1, 1, 0}));
    CHECK(isos == Reals({42.0}));
  }
  {  // triangle expected
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 2, LOs({0, 1, 2}), 3);
    mesh.add_tag(VERT, "coordinates", 2, OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals({0, 0, 1, -1, 1, 1}));
    auto isos = get_pad_isos(&mesh, 2, 1.0, 42.0, Read<I8>({1, 1, 0}));
    CHECK(are_close(isos, Reals({1.0})));
  }
  {  // tet with two bridges
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(VERT, "coordinates", 3, OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals(3 * 4, 0.0));
    auto isos = get_pad_isos(&mesh, 3, 1.0, 42.0, Read<I8>({1, 1, 0, 0, 0, 0}));
    CHECK(are_close(isos, Reals({42.0})));
  }
  {  // tet with three bridges, off-center
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(VERT, "coordinates", 3, OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals({0, 0, 0, 1, -1, 1, 1, 1, 1, 1, 0, 2}));
    auto isos = get_pad_isos(&mesh, 3, 1.0, 42.0, Read<I8>({1, 1, 1, 0, 0, 0}));
    CHECK(are_close(isos, Reals({42.0})));
  }
  {  // tet with three bridges, expected
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(VERT, "coordinates", 3, OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals({0, 0, 0, 1, -1, -1, 1, 1, -1, 1, 0, 2}));
    auto isos = get_pad_isos(&mesh, 3, 1.0, 42.0, Read<I8>({1, 1, 1, 0, 0, 0}));
    CHECK(are_close(isos, Reals({1.0})));
  }
  {  // edge-edge tet, off center
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(VERT, "coordinates", 3, OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals({0, 0, 0, 1, 0, 0, -1, 1, 0, -1, 1, 1}));
    auto isos = get_pad_isos(&mesh, 3, 1.0, 42.0, Read<I8>({0, 1, 1, 1, 1, 0}));
    CHECK(are_close(isos, Reals({42.0})));
  }
  {  // edge-edge tet, expected
    Mesh mesh(lib);
    build_from_elems2verts(&mesh, 3, LOs({0, 1, 2, 3}), 4);
    mesh.add_tag(VERT, "coordinates", 3, OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, Reals({0, 0, 0, 2, 0, 0, 1, 1, -1, 1, 1, 1}));
    auto isos = get_pad_isos(&mesh, 3, 1.0, 42.0, Read<I8>({0, 1, 1, 1, 1, 0}));
    CHECK(are_close(isos, Reals({1.0})));
  }
}

static void test_motion(Library* lib) {
  Mesh mesh(lib);
  build_box(&mesh, 1, 1, 0, 2, 2, 0);
  classify_by_angles(&mesh, Omega_h::PI / 4);
  auto sizes = find_implied_size(&mesh);
  mesh.add_tag(VERT, "size", 1, OMEGA_H_SIZE, OMEGA_H_DO_OUTPUT, sizes);
  auto coords_w = deep_copy(mesh.coords());
  coords_w.set(4 * 2 + 0, 0.74);
  coords_w.set(4 * 2 + 1, 0.26);
  mesh.set_coords(coords_w);
  CHECK(mesh.min_quality() < 0.5);
  AdaptOpts opts(&mesh);
  auto cands2verts = LOs({4});
  auto choices = get_motion_choices(&mesh, opts, cands2verts);
  CHECK(choices.cands_did_move.get(0));
  CHECK(choices.quals.get(0) > 0.85);
  auto verts_are_keys = Read<I8>({0, 0, 0, 0, 1, 0, 0, 0, 0});
  unpack_linearized_fields(&mesh, &mesh, choices.new_sol, verts_are_keys);
  CHECK(mesh.min_quality() == choices.quals.get(0));
}

static void test_find_last() {
  auto a = LOs({0, 3, 55, 12});
  CHECK(find_last(a, 98) < 0);
  CHECK(find_last(a, 12) == 3);
  CHECK(find_last(a, 55) == 2);
  CHECK(find_last(a, 3) == 1);
  CHECK(find_last(a, 0) == 0);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  test_edge_length();
  test_cubic();
  test_form_ortho_basis();
  test_qr_decomps();
  test_eigen_quadratic();
  test_eigen_cubic();
  test_least_squares();
  test_int128();
  test_repro_sum();
  test_sort();
  test_scan();
  test_intersect_metrics();
  test_fan_and_funnel();
  test_permute();
  test_invert_map();
  test_invert_adj();
  test_tri_align();
  test_form_uses();
  test_reflect_down();
  test_find_unique();
  test_hilbert();
  test_bbox();
  test_build_from_elems2verts(&lib);
  test_star(&lib);
  test_injective_map();
  test_dual(&lib);
  test_quality();
  test_file_components();
  test_linpart();
  test_expand();
  test_inertial_bisect(&lib);
  test_average_field(&lib);
  test_positivize();
  test_refine_qualities(&lib);
  test_mark_up_down(&lib);
  test_compare_meshes(&lib);
  test_swap2d_topology(&lib);
  test_swap3d_loop(&lib);
  test_file(&lib);
  test_xml();
  test_read_vtu(&lib);
  test_interpolate_isos();
  test_interpolate_metrics();
  test_element_implied_metric();
  test_recover_hessians(&lib);
  test_sf_scale(&lib);
  test_categorize_graph();
  test_circumcenter();
  test_lie();
  test_proximity(&lib);
  test_motion(&lib);
  test_find_last();
  CHECK(get_current_bytes() == 0);
}
