#include "internal.hpp"

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
  bool ok = decompose_eigen(m, q, l);
  CHECK(ok);
  CHECK(are_close(q,q_expect));
  CHECK(are_close(l,l_expect));
}

static void test_eigen_cubic(Matrix<3,3> m,
    Vector<3> l_expect) {
  Matrix<3,3> q;
  Vector<3> l;
  bool ok = decompose_eigen(m, q, l);
  CHECK(ok);
  CHECK(are_close(l,l_expect,1e-8,1e-8));
  CHECK(are_close(m, compose_eigen(q, l)));
}

static void test_eigen_cubic_ortho(Matrix<3,3> m,
    Vector<3> l_expect) {
  Matrix<3,3> q;
  Vector<3> l;
  bool ok = decompose_eigen(m, q, l);
  CHECK(ok);
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
  LOs scanned = offset_scan<LO>(LOs(3,1));
  CHECK(scanned == Read<LO>(4, 0, 1));
  }
  {
  LOs scanned = offset_scan<LO>(Read<I8>(3,1));
  CHECK(scanned == Read<LO>(4, 0, 1));
  }
}

static void test_invert_funnel() {
  CHECK(invert_funnel(LOs({0,0,1,1,2,2}), 3)
      == LOs({0,2,4,6}));
  CHECK(invert_funnel(LOs({0,0,0,2,2,2}), 3)
      == LOs({0,3,3,6}));
  CHECK(invert_funnel(LOs({0,0,0,0,0,0}), 3)
      == LOs({0,6,6,6}));
  CHECK(invert_funnel(LOs({2,2,2,2,2,2}), 3)
      == LOs({0,0,0,6}));
}

static void test_permute() {
  Reals data({0.1,0.2,0.3,0.4});
  LOs perm({3,2,1,0});
  Reals permuted = permute(perm, data);
  CHECK(permuted == Reals({0.4,0.3,0.2,0.1}));
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

static void test_invert_adj(map::InvertMethod method) {
  Adj tris2verts(LOs({0,1,2,2,3,0}));
  Read<GO> tri_globals({0,1});
  Adj verts2tris = invert(tris2verts, 3, 4, tri_globals, method);
  CHECK(verts2tris.a2ab == offset_scan<LO>(LOs({2,1,2,1})));
  CHECK(verts2tris.ab2b == LOs({0,1, 0, 0,1, 1}));
  CHECK(verts2tris.codes == Read<I8>({
        make_code(0, 0, 0),
        make_code(0, 0, 2),
        make_code(0, 0, 1),
        make_code(0, 0, 2),
        make_code(0, 0, 0),
        make_code(0, 0, 1)}));
}

static void test_invert_adj() {
  test_invert_adj(map::BY_SORTING);
  test_invert_adj(map::BY_ATOMICS);
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

static void test_reflect_down() {
  LOs eu2e;
  Read<I8> codes;
  {
  reflect_down<2>(LOs({}), LOs({}), eu2e, codes);
  CHECK(eu2e.size() == 0);
  CHECK(codes.size() == 0);
  }
  {
  reflect_down<2>(LOs({0,1}), LOs({0,1}), eu2e, codes);
  CHECK(eu2e == LOs({0}));
  CHECK(codes == Read<I8>({0}));
  }
  {
  reflect_down<2>(LOs({1,0}), LOs({0,1}), eu2e, codes);
  CHECK(eu2e == LOs({0}));
  CHECK(codes == Read<I8>({make_code(false, 1, 0)}));
  }
  {
  reflect_down<2>(LOs({1,0,0,1}), LOs({0,1}), eu2e, codes);
  CHECK(eu2e == LOs({0,0}));
  CHECK(codes == Read<I8>({
        make_code(false, 1, 0),
        make_code(false, 0, 0)}));
  }
  {
  reflect_down<3>(LOs({0,1,2,0,2,1}), LOs({0,1,2}), eu2e, codes);
  CHECK(eu2e == LOs({0,0}));
  CHECK(codes == Read<I8>({
        make_code(false, 0, 0),
        make_code(true, 0, 0)}));
  }
  {
  reflect_down<3>(LOs({0,1,2,0,2,1,1,2,0}), LOs({0,1,2}), eu2e, codes);
  CHECK(eu2e == LOs({0,0,0}));
  CHECK(codes == Read<I8>({
        make_code(false, 0, 0),
        make_code(true, 0, 0),
        make_code(false, 2, 0)}));
  }
  {
  reflect_down<3>(LOs({0,1,2,2,3,0}), LOs({0,1,2,2,3,0}), eu2e, codes);
  CHECK(eu2e == LOs({0,1}));
  CHECK(codes == Read<I8>({0,0}));
  }
  {
  reflect_down<3>(LOs({0,1,2, 2,0,3, 2,3,0}), LOs({2,3,0,0,1,2}), eu2e, codes);
  CHECK(eu2e == LOs({1,0,0}));
  CHECK(codes == Read<I8>({0,make_code(true,0,0),0}));
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

int main(int argc, char** argv) {
  init(argc, argv);
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
  test_invert_funnel();
  test_permute();
  test_invert_map();
  test_invert_adj();
  test_tri_align();
  test_reflect_down();
  test_form_uses();
  fini();
}
