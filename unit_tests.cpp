#include "internal.hpp"

using namespace osh;

template <Int m, Int n>
static void test_qr_decomp(Matrix<m,n> a) {
  Matrix<m,n> q;
  Matrix<n,n> r;
  decompose_qr_reduced(a, q, r);
  CHECK(are_close(a, q * r));
  CHECK(are_close(transpose(q) * q, identity_matrix<n,n>()));
}

static void test_qr_decomps() {
  test_qr_decomp(Matrix<3,3>({
       0, 0, 0,
       0, 0, 0,
       0, 0, 0}));
  test_qr_decomp(Matrix<3,3>({
       EPSILON, 0, 0,
       EPSILON, EPSILON, 0,
       EPSILON, EPSILON, EPSILON}));
  test_qr_decomp(Matrix<3,3>({
      12, -51,  4,
       6, 167,-68,
      -4,  24,-41}));
}

static void test_form_ortho_basis() {
  auto n = normalize(vector_3(1,1,1));
  auto f = form_ortho_basis(n);
  CHECK(are_close(f[0], n));
  CHECK(are_close(transpose(f) * f, identity_matrix<3,3>()));
}

static void test_least_squares() {
  Matrix<4,2> m({
      1, 1,
      1, 2,
      1, 3,
      1, 4});
  Vector<4> b({6, 5, 7, 10});
  Vector<2> x;
  CHECK(solve_least_squares_qr(m, b, x));
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
  Reals a({std::exp2(int(20)),std::exp2(int(-20))});
  Real sum = repro_sum(a);
  CHECK(sum == std::exp2(20) + std::exp2(int(-20)));
}

static void test_cubic(Real a, Real b, Real c,
    Int nroots_wanted, Few<Real, 3> roots_wanted,
    Few<Int, 3> mults_wanted) {
  Few<Real, 3> roots;
  Few<Int, 3> mults;
  Int nroots = solve_cubic(a, b, c, roots, mults);
  CHECK(nroots == nroots_wanted);
  for (Int i = 0; i < nroots; ++i) {
    CHECK(mults[i] == mults_wanted[i]);
    CHECK(are_close(roots[i], roots_wanted[i]));
  }
}

static void test_cubic() {
  test_cubic(0, 0, 0,
      1, Few<Real,3>({0}), Few<Int,3>({3}));
  test_cubic(-3. / 2., -3. / 2., 1.,
      3, Few<Real,3>({2,-1,.5}), Few<Int,3>({1,1,1}));
  test_cubic(0, -3., 2.,
      2, Few<Real,3>({-2,1}), Few<Int,3>({1,2}));
  test_cubic(3, -6, -8,
      3, Few<Real,3>({2,-4,-1}), Few<Int,3>({1,1,1}));
}

static void test_eigen_cubic(Matrix<3,3> m,
    Matrix<3,3> q_expect, Vector<3> l_expect) {
  Matrix<3,3> q;
  Vector<3> l;
  decompose_eigen(m, q, l);
  CHECK(are_close(q,q_expect));
  CHECK(are_close(l,l_expect));
}

static void test_eigen_cubic(Matrix<3,3> m,
    Vector<3> l_expect) {
  Matrix<3,3> q;
  Vector<3> l;
  decompose_eigen(m, q, l);
  CHECK(are_close(l,l_expect,1e-8,1e-8));
  CHECK(are_close(m, compose_eigen(q, l)));
}

static void test_eigen_cubic_ortho(Matrix<3,3> m,
    Vector<3> l_expect) {
  Matrix<3,3> q;
  Vector<3> l;
  decompose_eigen(m, q, l);
  CHECK(are_close(transpose(q) * q, identity_matrix<3,3>(), 1e-8, 1e-8));
  CHECK(are_close(l,l_expect, 1e-8, 1e-8));
  CHECK(are_close(m, compose_ortho(q, l), 1e-8, 1e-8));
}

static void test_eigen_metric(Vector<3> h) {
  auto q = rotate(PI / 4., vector_3(0,0,1)) *
           rotate(PI / 4., vector_3(0,1,0));
  CHECK(are_close(transpose(q) * q, identity_matrix<3,3>()));
  auto l = metric_eigenvalues(h);
  auto a = compose_ortho(q, l);
  test_eigen_cubic_ortho(a, l);
}

static void test_eigen_cubic() {
  test_eigen_cubic(
      identity_matrix<3,3>(),
      identity_matrix<3,3>(),
      vector_3(1,1,1));
  test_eigen_cubic(
      matrix_3x3(0,0,0,0,0,0,0,0,0),
      identity_matrix<3,3>(),
      vector_3(0,0,0));
  test_eigen_cubic(
      matrix_3x3(
        -1, 3, -1,
        -3, 5, -1,
        -3, 3,  1),
      vector_3(1,2,2));
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
    Vector<3> h1,
    Vector<3> h2,
    Vector<3> hi_expect) {
  auto q = rotate(PI / 4., vector_3(0,0,1)) *
           rotate(PI / 4., vector_3(0,1,0));
  auto m1 = compose_metric(q, h1);
  auto m2 = compose_metric(q, h2);
  auto mi = intersect_metrics(m1, m2);
  /* if we decompose it, the eigenvectors may
     get re-ordered. */
  for (Int i = 0; i < 3; ++i) {
    CHECK(are_close(metric_desired_length(mi, q[i]),
                    hi_expect[i], 1e-3));
  }
}

static void test_intersect_metrics() {
  test_intersect_ortho_metrics(
      vector_3(0.5,   1, 1),
      vector_3(  1, 0.5, 1),
      vector_3(0.5, 0.5, 1));
  test_intersect_ortho_metrics(
      vector_3(1e-3, 1,    1),
      vector_3(   1, 1, 1e-3),
      vector_3(1e-3, 1, 1e-3));
  test_intersect_ortho_metrics(
      vector_3(1e-3, 1e-3,    1),
      vector_3(   1,    1, 1e-3),
      vector_3(1e-3, 1e-3, 1e-3));
  test_intersect_ortho_metrics(
      vector_3(1e-6, 1e-3, 1e-3),
      vector_3(1e-3, 1e-3, 1e-6),
      vector_3(1e-6, 1e-3, 1e-6));
}

static void test_sort() {
  {
  LOs a({0,1});
  LOs perm = sort_by_keys(a);
  CHECK(perm == LOs({0,1}));
  }
  {
  LOs a({0,2,0,1});
  LOs perm = sort_by_keys<LO,2>(a);
  CHECK(perm == LOs({1,0}));
  }
  {
  LOs a({0,2,1,1});
  LOs perm = sort_by_keys<LO,2>(a);
  CHECK(perm == LOs({0,1}));
  }
  {
  LOs a({1,2,3,1,2,2,3,0,0});
  LOs perm = sort_by_keys<LO,3>(a);
  CHECK(perm == LOs({1,0,2}));
  }
}

static void test_scan() {
  {
  LOs scanned = offset_scan(LOs(3,1));
  CHECK(scanned == Read<LO>(4, 0, 1));
  }
  {
  LOs scanned = offset_scan(Read<I8>(3,1));
  CHECK(scanned == Read<LO>(4, 0, 1));
  }
}

static void test_fan_and_funnel() {
  CHECK(invert_funnel(LOs({0,0,1,1,2,2}), 3)
      == LOs({0,2,4,6}));
  CHECK(invert_fan(LOs({0,2,4,6}))
      == LOs({0,0,1,1,2,2}));
  CHECK(invert_funnel(LOs({0,0,0,2,2,2}), 3)
      == LOs({0,3,3,6}));
  CHECK(invert_fan(LOs({0,3,3,6}))
      == LOs({0,0,0,2,2,2}));
  CHECK(invert_funnel(LOs({0,0,0,0,0,0}), 3)
      == LOs({0,6,6,6}));
  CHECK(invert_fan(LOs({0,6,6,6}))
      == LOs({0,0,0,0,0,0}));
  CHECK(invert_funnel(LOs({2,2,2,2,2,2}), 3)
      == LOs({0,0,0,6}));
  CHECK(invert_fan(LOs({0,0,0,6}))
      == LOs({2,2,2,2,2,2}));
}

static void test_permute() {
  Reals data({0.1,0.2,0.3,0.4});
  LOs perm({3,2,1,0});
  Reals permuted = unmap(perm, data, 1);
  CHECK(permuted == Reals({0.4,0.3,0.2,0.1}));
  Reals back = permute(permuted, perm, 1);
  CHECK(back == data);
}

// these tests can have degree at most 1
// because map::invert doesn't have to be
// deterministic in local ordering
static void test_invert_map(map::InvertMethod method) {
  {
  LOs hl2l({});
  LOs l2lh;
  LOs lh2hl;
  map::invert(hl2l, 4, l2lh, lh2hl, method);
  CHECK(l2lh == LOs(5,0));
  CHECK(lh2hl == LOs({}));
  }
  {
  LOs hl2l({0,1,2,3});
  LOs l2lh;
  LOs lh2hl;
  map::invert(hl2l, 4, l2lh, lh2hl, method);
  CHECK(l2lh == LOs(5,0,1));
  CHECK(lh2hl == LOs(4,0,1));
  }
}

static void test_invert_map() {
  test_invert_map(map::BY_SORTING);
  test_invert_map(map::BY_ATOMICS);
}

static void test_invert_adj() {
  Adj tris2verts(LOs({0,1,2,2,3,0}));
  Read<GO> tri_globals({0,1});
  Adj verts2tris = invert(tris2verts, 3, 4, tri_globals);
  CHECK(verts2tris.a2ab == offset_scan(LOs({2,1,2,1})));
  CHECK(verts2tris.ab2b == LOs({0,1, 0, 0,1, 1}));
  CHECK(verts2tris.codes == Read<I8>({
        make_code(0, 0, 0),
        make_code(0, 0, 2),
        make_code(0, 0, 1),
        make_code(0, 0, 2),
        make_code(0, 0, 0),
        make_code(0, 0, 1)}));
}

static bool same_adj(Int a[], Int b[]) {
  for (Int i = 0; i < 3; ++i)
    if (a[i] != b[i])
      return false;
  return true;
}

static void test_tri_align() {
  Int ident[3] = {0,1,2};
  Int out[3];
  Int out2[3];
  /* check that flipping and rotating do what we want */
  {
  align_adj<3,Int>(make_code(true, 0, 0), ident, out);
  Int expect[3] = {0,2,1};
  CHECK(same_adj(out, expect));
  }
  {
  align_adj<3>(make_code(false, 1, 0), ident, out);
  Int expect[3] = {2,0,1};
  CHECK(same_adj(out, expect));
  }
  {
  align_adj<3>(make_code(false, 2, 0), ident, out);
  Int expect[3] = {1,2,0};
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
  CHECK(form_uses(LOs({0,1,2}),2,1) ==
      LOs({0,1,1,2,2,0}));
  CHECK(form_uses(LOs({0,1,2,3}),3,1) ==
      LOs({0,1,1,2,2,0,0,3,1,3,2,3}));
  CHECK(form_uses(LOs({0,1,2,3}),3,2) ==
      LOs({0,2,1,0,1,3,1,2,3,2,0,3}));
}

static void test_reflect_down() {
  Adj a;
  a = reflect_down(LOs({}),LOs({}),0,2,1);
  CHECK(a.ab2b == LOs({}));
  CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({}),LOs({}),0,3,1);
  CHECK(a.ab2b == LOs({}));
  CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({}),LOs({}),0,3,2);
  CHECK(a.ab2b == LOs({}));
  CHECK(a.codes == Read<I8>({}));
  a = reflect_down(LOs({0,1,2}),LOs({0,1,1,2,2,0}),3,2,1);
  CHECK(a.ab2b == LOs({0,1,2}));
  CHECK(a.codes == Read<I8>({0,0,0}));
  a = reflect_down(LOs({0,1,2,3}),
      LOs({0,1,1,2,2,0,0,3,1,3,2,3}),4,3,1);
  CHECK(a.ab2b == LOs({0,1,2,3,4,5}));
  CHECK(a.codes == Read<I8>({0,0,0,0,0,0}));
  a = reflect_down(LOs({0,1,2,3}),
      LOs({0,2,1,0,1,3,1,2,3,2,0,3}),4,3,2);
  CHECK(a.ab2b == LOs({0,1,2,3}));
  CHECK(a.codes == Read<I8>({0,0,0,0}));
  a = reflect_down(LOs({0,1,2,3}),
      LOs({0,1,2,0,3,1,1,3,2,2,3,0}),4,3,2);
  CHECK(a.ab2b == LOs({0,1,2,3}));
  CHECK(a.codes == Read<I8>(4, make_code(true,0,0)));
  a = reflect_down(LOs({0,1,2,2,3,0}),
      LOs({0,1,1,2,2,3,3,0,0,2}),4,2,1);
  CHECK(a.ab2b == LOs({0,1,4,2,3,4}));
}

static void test_find_unique() {
  CHECK(find_unique(LOs({}), 2, 1) == LOs({}));
  CHECK(find_unique(LOs({}), 3, 1) == LOs({}));
  CHECK(find_unique(LOs({}), 3, 2) == LOs({}));
  CHECK(find_unique(LOs({0,1,2,2,3,0}),2,1) ==
      LOs({0,1,0,2,3,0,1,2,2,3}));
}

static void test_hilbert() {
  /* this is the original test from Skilling's paper */
  hilbert::coord_t X[3] = {5,10,20}; // any position in 32x32x32 cube
  hilbert::AxestoTranspose(X, 5, 3); // Hilbert transpose for 5 bits and 3 dimensions
  std::stringstream stream;
  stream << "Hilbert integer = "
    << (X[0]>>4 & 1) << (X[1]>>4 & 1) << (X[2]>>4 & 1)
    << (X[0]>>3 & 1) << (X[1]>>3 & 1) << (X[2]>>3 & 1)
    << (X[0]>>2 & 1) << (X[1]>>2 & 1) << (X[2]>>2 & 1)
    << (X[0]>>1 & 1) << (X[1]>>1 & 1) << (X[2]>>1 & 1)
    << (X[0]>>0 & 1) << (X[1]>>0 & 1) << (X[2]>>0 & 1)
    << " = 7865 check";
  std::string expected =
  "Hilbert integer = 001111010111001 = 7865 check";
  CHECK(stream.str() == expected);
  hilbert::coord_t Y[3];
  hilbert::untranspose(X, Y, 5, 3);
  std::stringstream stream2;
  stream2 << "Hilbert integer = "
    <<(Y[0]>>4 & 1)<<(Y[0]>>3 & 1)<<(Y[0]>>2 & 1)<<(Y[0]>>1 & 1)<<(Y[0]>>0 & 1)
    <<(Y[1]>>4 & 1)<<(Y[1]>>3 & 1)<<(Y[1]>>2 & 1)<<(Y[1]>>1 & 1)<<(Y[1]>>0 & 1)
    <<(Y[2]>>4 & 1)<<(Y[2]>>3 & 1)<<(Y[2]>>2 & 1)<<(Y[2]>>1 & 1)<<(Y[2]>>0 & 1)
    << " = 7865 check";
  CHECK(stream2.str() == expected);
}

static void test_bbox() {
  CHECK(are_close(BBox<2>(vector_2(-3,-3),vector_2(3,3)),
        find_bounding_box<2>(Reals({0,-3,3,0,0,3,-3,0}))));
  CHECK(are_close(BBox<3>(vector_3(-3,-3,-3),vector_3(3,3,3)),
        find_bounding_box<3>(Reals({
            0,-3, 0,
            3, 0, 0,
            0, 3, 0,
           -3, 0, 0,
            0, 0,-3,
            0, 0, 3}))));
}

static void test_build_from_elems2verts() {
  {
  Mesh mesh;
  build_from_elems2verts(mesh,2,LOs({0,1,2}),3);
  CHECK(mesh.ask_down(2,0).ab2b == LOs({0,1,2}));
  CHECK(mesh.ask_down(2,1).ab2b == LOs({0,2,1}));
  CHECK(mesh.ask_down(1,0).ab2b ==
      LOs({0,1,2,0,1,2}));
  }
  {
  Mesh mesh;
  build_from_elems2verts(mesh,3,LOs({0,1,2,3}),4);
  CHECK(mesh.ask_down(3,0).ab2b == LOs({0,1,2,3}));
  }
}

static void test_star() {
  {
  Mesh mesh;
  build_from_elems2verts(mesh,2,LOs({0,1,2}),3);
  Adj v2v = mesh.ask_star(VERT);
  CHECK(v2v.a2ab == LOs(4,0,2));
  CHECK(v2v.ab2b == LOs({1,2,0,2,0,1}));
  Adj e2e = mesh.ask_star(EDGE);
  CHECK(e2e.a2ab == LOs(4,0,2));
  CHECK(e2e.ab2b == LOs({2,1,0,2,1,0}));
  }
  {
  Mesh mesh;
  build_from_elems2verts(mesh,3,LOs({0,1,2,3}),4);
  Adj v2v = mesh.ask_star(VERT);
  CHECK(v2v.a2ab == LOs(5,0,3));
  CHECK(v2v.ab2b == LOs({
        1,2,3,0,2,3,0,1,3,0,1,2}));
  Adj e2e = mesh.ask_star(EDGE);
  CHECK(e2e.a2ab == LOs(7,0,5));
  CHECK(e2e.ab2b == LOs({
        1,3,4,2,5,
        3,0,2,5,4,
        0,4,5,1,3,
        0,1,5,4,2,
        2,0,3,5,1,
        1,2,4,3,0}));
  }
}

static void test_injective_map() {
  LOs primes2ints({2,3,5,7});
  LOs ints2primes = invert_injective_map(primes2ints, 8);
  CHECK(ints2primes == LOs({-1,-1,0,1,-1,2,-1,3}));
}

static void test_dual() {
  Mesh mesh;
  build_from_elems2verts(mesh,2,LOs({0,1,2,2,3,0}),4);
  auto t2t = mesh.ask_dual();
  auto t2tt = t2t.a2ab;
  auto tt2t = t2t.ab2b;
  CHECK(t2tt == offset_scan(LOs({1,1})));
  CHECK(tt2t == LOs({1,0}));
}

static void test_quality() {
  Few<Vector<2>, 3> perfect_tri({
      vector_2( 1, 0),
      vector_2( 0, sqrt(3.0)),
      vector_2(-1, 0)});
  Few<Vector<3>, 4> perfect_tet({
      vector_3( 1, 0, -1.0 / sqrt(2.0)),
      vector_3(-1, 0, -1.0 / sqrt(2.0)),
      vector_3( 0,-1,  1.0 / sqrt(2.0)),
      vector_3( 0, 1,  1.0 / sqrt(2.0))});
  Few<Vector<2>, 3> flat_tri({
      vector_2( 1, 0),
      vector_2( 0, 0),
      vector_2(-1, 0)});
  Few<Vector<3>, 4> flat_tet({
      vector_3( 1, 0, 0),
      vector_3(-1, 0, 0),
      vector_3( 0,-1, 0),
      vector_3( 0, 1, 0)});
  Few<Vector<2>, 3> inv_tri({
      vector_2( 1, 0),
      vector_2(-1, 0),
      vector_2( 0, sqrt(3.0))});
  Few<Vector<3>, 4> inv_tet({
      vector_3( 1, 0, -1.0 / sqrt(2.0)),
      vector_3(-1, 0, -1.0 / sqrt(2.0)),
      vector_3( 0, 1,  1.0 / sqrt(2.0)),
      vector_3( 0,-1,  1.0 / sqrt(2.0))});
  Matrix<2,2> id_metric_2 = identity_matrix<2,2>();
  Matrix<3,3> id_metric_3 = identity_matrix<3,3>();
  Matrix<2,2> x_metric_2 = compose_metric(identity_matrix<2,2>(),
      vector_2(1, 0.5));
  Matrix<3,3> x_metric_3 = compose_metric(identity_matrix<3,3>(),
      vector_3(1, 1, 0.5));
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
  CHECK(are_close(real_triangle_quality(perfect_tri), 1.0));
  CHECK(are_close(real_tet_quality(perfect_tet), 1.0));
  CHECK(are_close(real_triangle_quality(flat_tri), 0.0));
  CHECK(are_close(real_tet_quality(flat_tet), 0.0));
  CHECK(real_triangle_quality(inv_tri) < 0.0);
  CHECK(real_tet_quality(inv_tet) < 0.0);
  CHECK(are_close(metric_triangle_quality(perfect_tri, id_metric_2), 1.0));
  CHECK(are_close(metric_tet_quality(perfect_tet, id_metric_3), 1.0));
  CHECK(are_close(metric_triangle_quality(flat_tri, id_metric_2), 0.0));
  CHECK(are_close(metric_tet_quality(flat_tet, id_metric_3), 0.0));
  CHECK(metric_triangle_quality(inv_tri, id_metric_2) < 0.0);
  CHECK(metric_tet_quality(inv_tet, id_metric_3) < 0.0);
  CHECK(are_close(metric_triangle_quality(x_tri, x_metric_2), 1.0));
  CHECK(are_close(metric_tet_quality(x_tet, x_metric_3), 1.0));
}

static void test_file() {
  using namespace binary;
  std::stringstream stream;
  std::string s = "foo";
  LO n = 10;
#ifdef OSH_USE_ZLIB
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
  Read<GO> globals({6,5,4,3,2,1,0});
  auto remotes = globals_to_linear_owners(globals, total, comm_size);
  CHECK(remotes.ranks == Read<I32>({1,1,1,0,0,0,0}));
  CHECK(remotes.idxs == Read<I32>({2,1,0,3,2,1,0}));
}

static void test_expand() {
  auto fan = offset_scan(LOs({2,1,3}));
  Reals data({2.2,3.14,42.0});
  CHECK(expand(data,fan,1) ==
      Reals({2.2,2.2,3.14,42.0,42.0,42.0}));
}

static void test_inertial_bisect() {
  Reals coords({ 2, 1,0,
                 2,-1,0,
                -2, 1,0,
                -2,-1,0});
  Reals masses(4, 1);
  auto self = Comm::self();
  Real tolerance = 0.0;
  Vector<3> axis;
  auto marked = inertia::mark_bisection(self, coords, masses, tolerance,
      axis);
  CHECK(marked == Read<I8>({1,1,0,0}));
  marked = inertia::mark_bisection_given_axis(self, coords, masses, tolerance,
      vector_3(0,1,0));
  CHECK(marked == Read<I8>({1,0,1,0}));
}

static void test_average_field() {
  Mesh mesh;
  build_box(mesh, 1, 1, 0, 1, 1, 0);
  Reals v2x({2,1,3,2});
  auto e2x = average_field(mesh, 2, LOs({0,1}), 1, v2x);
  CHECK(are_close(e2x, Reals({5. / 3., 7. / 3.})));
}

template <Int n>
static void test_positivize(Vector<n> pos) {
  auto neg = pos * -1.0;
  CHECK(are_close(positivize(pos), pos));
  CHECK(are_close(positivize(neg), pos));
}

static void test_positivize() {
  test_positivize(vector_3( 1, 1, 1));
  test_positivize(vector_3( 1,-1, 1));
  test_positivize(vector_2(-1, 1));
  test_positivize(vector_2( 1, 1));
}

static void test_refine_qualities() {
  Mesh mesh;
  build_box(mesh, 1, 1, 0, 1, 1, 0);
  LOs candidates(mesh.nents(EDGE), 0, 1);
  auto quals = refine_qualities(mesh, candidates);
  CHECK(are_close(quals,
        Reals({0.494872,0.494872,0.866025,0.494872,0.494872}),
        1e-4));
  mesh.add_tag(VERT, "metric", symm_dofs(2), OSH_METRIC,
      repeat_symm(mesh.nverts(), identity_matrix<2,2>()));
  auto quals2 = refine_qualities(mesh, candidates);
  CHECK(are_close(quals2, quals));
}

static void test_i8_array_print() {
  Read<I8> a({1,1,1});
  std::stringstream stream;
  stream << a;
  std::string str = stream.str();
  CHECK(str == " 1 1 1\n");
}

static void test_mark_up_down() {
  Mesh mesh;
  build_box(mesh, 1, 1, 0, 1, 1, 0);
  CHECK(mark_down(mesh, TRI, VERT, Read<I8>({1,0}))
        == Read<I8>({1,1,0,1}));
  CHECK(mark_up(mesh, VERT, TRI, Read<I8>({0,1,0,0}))
        == Read<I8>({1,0}));
}

static void test_compare_meshes() {
  Mesh a;
  build_box(a, 1, 1, 0, 4, 4, 0);
  CHECK(a == a);
  Mesh b = a;
  b.reorder();
  CHECK(a == b);
  b.add_tag<I8>(VERT, "foo", 1, OSH_DONT_TRANSFER,
      Read<I8>(b.nverts(), 1));
  CHECK(!(a == b));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  test_cubic();
  test_form_ortho_basis();
  test_qr_decomps();
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
  test_build_from_elems2verts();
  test_star();
  test_injective_map();
  test_dual();
  test_quality();
  test_file();
  test_linpart();
  test_expand();
  test_inertial_bisect();
  test_average_field();
  test_positivize();
  test_refine_qualities();
  test_i8_array_print();
  test_mark_up_down();
  test_compare_meshes();
}
