#include "internal.hpp"

template <UInt m, UInt n>
static void test_qr_decomp(Matrix<m,n> a) {
  Matrix<m,n> q;
  Matrix<n,n> r;
  decompose_qr_reduced(a, q, r);
  CHECK(are_close(a, q * r));
  CHECK(are_close(transpose(q) * q, identity_matrix<m,n>()));
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

static void test_eigen_decomp() {
  auto q = rotate(PI / 4, vector_3(0,0,1)) *
           rotate(PI / 4, vector_3(0,1,0));
  CHECK(are_close(transpose(q) * q, identity_matrix<3,3>()));
  auto l = matrix_3x3(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1e-6);
  auto a = q * l * transpose(q);
  Matrix<3,3> q2;
  Matrix<3,3> l2;
  decompose_eigen_qr(a, q2, l2);
  CHECK(are_close(transpose(q2) * q2, identity_matrix<3,3>()));
  CHECK(are_close(q2 * l2 * transpose(q2), a));
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
    UInt nroots_wanted, Few<Real, 3> roots_wanted,
    Few<UInt, 3> mults_wanted) {
  std::cerr << "test_cubic\n";
  Few<Real, 3> roots;
  Few<UInt, 3> mults;
  UInt nroots = solve_cubic(a, b, c, &roots[0], &mults[0]);
  std::cerr << "nroots " << nroots << " wanted " << nroots_wanted << '\n';
  CHECK(nroots == nroots_wanted);
  for (UInt i = 0; i < nroots; ++i) {
    std::cerr << "root " << roots[i] << " wanted " << roots_wanted[i] << '\n';
    CHECK(mults[i] == mults_wanted[i]);
    CHECK(are_close(roots[i], roots_wanted[i]));
  }
}

static void test_cubic() {
  Few<Real, 3> roots;
  Few<UInt, 3> mults;
  roots[0] = 0;
  mults[0] = 3;
  test_cubic(0, 0, 0,
      1, roots, mults);
  roots[0] = -1;
  roots[1] =  2;
  roots[2] = 0.5;
  mults[0] = 1;
  mults[1] = 1;
  mults[2] = 1;
  test_cubic(-3. / 2., -3. / 2., 1.,
      3, roots, mults);
  roots[0] = -2;
  roots[1] =  1;
  mults[0] = 1;
  mults[1] = 2;
  test_cubic(0, -3., 2.,
      2, roots, mults);
}

int main(int argc, char** argv) {
  init(argc, argv);
  test_form_ortho_basis();
  test_qr_decomps();
  test_eigen_decomp();
  test_least_squares();
  test_int128();
  test_repro_sum();
  test_cubic();
  fini();
}
