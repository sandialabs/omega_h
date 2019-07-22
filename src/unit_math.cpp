#include <Omega_h_array_ops.hpp>
#include <Omega_h_eigen.hpp>
#include <Omega_h_lie.hpp>
#include <Omega_h_metric_intersect.hpp>
#include <Omega_h_most_normal.hpp>
#include <Omega_h_shape.hpp>
#include <Omega_h_svd.hpp>

using namespace Omega_h;

static void test_edge_length() {
  OMEGA_H_CHECK(are_close(1., anisotropic_edge_length(1., 1.)));
  OMEGA_H_CHECK(anisotropic_edge_length(1., 2.) > 1.);
  OMEGA_H_CHECK(anisotropic_edge_length(1., 2.) < 1.5);
}

static void test_least_squares() {
  Matrix<4, 2> m({1, 1, 1, 2, 1, 3, 1, 4});
  Vector<4> b({6, 5, 7, 10});
  auto x = solve_using_qr(m, b);
  OMEGA_H_CHECK(are_close(x, vector_2(3.5, 1.4)));
}

static void test_power() {
  auto x = 3.14159;
  OMEGA_H_CHECK(x == power(x, 1, 1));
  OMEGA_H_CHECK(x == power(x, 2, 2));
  OMEGA_H_CHECK(x == power(x, 3, 3));
  OMEGA_H_CHECK(are_close(x * x, power(x, 2, 1)));
  OMEGA_H_CHECK(are_close(x * x * x, power(x, 3, 1)));
  OMEGA_H_CHECK(are_close(std::sqrt(x), power(x, 1, 2)));
  OMEGA_H_CHECK(are_close(std::cbrt(x), power(x, 1, 3)));
  OMEGA_H_CHECK(are_close(std::sqrt(x * x * x), power(x, 3, 2)));
  OMEGA_H_CHECK(are_close(std::cbrt(x * x), power(x, 2, 3)));
}

static void test_cubic(Few<Real, 3> coeffs, Int nroots_wanted,
    Few<Real, 3> roots_wanted, Few<Int, 3> mults_wanted) {
  auto roots = find_polynomial_roots(coeffs);
  OMEGA_H_CHECK(roots.n == nroots_wanted);
  for (Int i = 0; i < roots.n; ++i) {
    OMEGA_H_CHECK(roots.mults[i] == mults_wanted[i]);
    OMEGA_H_CHECK(are_close(roots.values[i], roots_wanted[i]));
  }
}

static void test_cubic() {
  test_cubic({0, 0, 0}, 1, Few<Real, 3>({0}), Few<Int, 3>({3}));
  test_cubic({1., -3. / 2., -3. / 2.}, 3, Few<Real, 3>({2, -1, .5}),
      Few<Int, 3>({1, 1, 1}));
  test_cubic({2., -3., 0}, 2, Few<Real, 3>({-2, 1}), Few<Int, 3>({1, 2}));
  test_cubic({-8, -6, 3}, 3, Few<Real, 3>({2, -4, -1}), Few<Int, 3>({1, 1, 1}));
}

static void test_form_ortho_basis() {
  auto n = normalize(vector_3(1, 1, 1));
  auto f = form_ortho_basis(n);
  OMEGA_H_CHECK(are_close(f[0], n));
  OMEGA_H_CHECK(are_close(transpose(f) * f, identity_matrix<3, 3>()));
}

template <Int m, Int n>
static void test_qr_decomp(Matrix<m, n> a) {
  auto qr = factorize_qr_householder(m, n, a);
  auto r = qr.r;
  auto q = identity_matrix<m, n>();
  for (Int j = 0; j < n; ++j) implicit_q_x(m, n, q[j], qr.v);
  OMEGA_H_CHECK(are_close(a, q * r));
  OMEGA_H_CHECK(are_close(transpose(q) * q, identity_matrix<n, n>()));
}

static void test_qr_decomps() {
  test_qr_decomp(identity_matrix<3, 3>());
  test_qr_decomp(Matrix<3, 3>({EPSILON, 0, 0, 0, EPSILON, 0, 0, 0, EPSILON}));
  test_qr_decomp(Matrix<3, 3>({12, -51, 4, 6, 167, -68, -4, 24, -41}));
}

template <Int dim>
static void test_eigen(
    Matrix<dim, dim> m, Matrix<dim, dim> q_expect, Vector<dim> l_expect) {
  auto ed = decompose_eigen(m);
  auto q = ed.q;
  auto l = ed.l;
  OMEGA_H_CHECK(are_close(q, q_expect));
  OMEGA_H_CHECK(are_close(l, l_expect));
}

template <Int dim>
static void test_eigen(Matrix<dim, dim> m, Vector<dim> l_expect) {
  auto ed = decompose_eigen(m);
  auto q = ed.q;
  auto l = ed.l;
  OMEGA_H_CHECK(are_close(l, l_expect, 1e-8, 1e-8));
  OMEGA_H_CHECK(are_close(m, compose_eigen(q, l)));
}

static void test_eigen_cubic_ortho(Matrix<3, 3> m, Vector<3> l_expect) {
  auto ed = decompose_eigen(m);
  auto q = ed.q;
  auto l = ed.l;
  OMEGA_H_CHECK(
      are_close(transpose(q) * q, identity_matrix<3, 3>(), 1e-8, 1e-8));
  OMEGA_H_CHECK(are_close(l, l_expect, 1e-8, 1e-8));
  OMEGA_H_CHECK(are_close(m, compose_ortho(q, l), 1e-8, 1e-8));
}

static void test_eigen_metric(Vector<3> h) {
  auto q =
      rotate(PI / 4., vector_3(0, 0, 1)) * rotate(PI / 4., vector_3(0, 1, 0));
  OMEGA_H_CHECK(are_close(transpose(q) * q, identity_matrix<3, 3>()));
  auto l = metric_eigenvalues_from_lengths(h);
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

template <Int dim>
static void test_eigen_jacobi(
    Matrix<dim, dim> a, Matrix<dim, dim> expect_q, Vector<dim> expect_l) {
  auto ed = decompose_eigen_jacobi(a);
  ed = sort_by_magnitude(ed);
  OMEGA_H_CHECK(are_close(ed.q, expect_q));
  OMEGA_H_CHECK(are_close(ed.l, expect_l));
}

static void test_eigen_jacobi_sign_bug() {
  auto sign_bug_input = matrix_3x3(0.99999999998511147, 5.3809065327405379e-11,
      9.7934015130085018e-10, 5.3809065327405379e-11, 0.99999999995912181,
      -1.676999986436999e-09, 9.7934015130085018e-10, -1.676999986436999e-09,
      0.99999995816580101);
  decompose_eigen_jacobi(sign_bug_input);
}

static void test_eigen_jacobi() {
  test_eigen_jacobi(
      identity_matrix<2, 2>(), identity_matrix<2, 2>(), vector_2(1, 1));
  test_eigen_jacobi(
      identity_matrix<3, 3>(), identity_matrix<3, 3>(), vector_3(1, 1, 1));
  test_eigen_jacobi(matrix_2x2(2, 1, 1, 2),
      matrix_2x2(1, 1, 1, -1) / std::sqrt(2), vector_2(3, 1));
  test_eigen_jacobi(matrix_3x3(2, 0, 0, 0, 3, 4, 0, 4, 9),
      Matrix<3, 3>({normalize(vector_3(0, 1, 2)), normalize(vector_3(1, 0, 0)),
          normalize(vector_3(0, 2, -1))}),
      vector_3(11, 2, 1));
  test_eigen_jacobi_sign_bug();
}

static void test_most_normal() {
  {
    Few<Vector<3>, 3> N = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    auto N_c = get_most_normal_normal(N, 3);
    OMEGA_H_CHECK(are_close(N_c, normalize(vector_3(1, 1, 1))));
  }
  {
    Few<Vector<3>, 4> N = {
        {1, 0, 0}, {0, 0, 1}, normalize(vector_3(1, 1, 1)), {0, 1, 0}};
    auto N_c = get_most_normal_normal(N, 4);
    OMEGA_H_CHECK(are_close(N_c, normalize(vector_3(1, 1, 1))));
  }
  {
    Few<Vector<3>, 3> N = {{1, 0, 0}, {0, 1, 0}, {0, 1, 0}};
    auto N_c = get_most_normal_normal(N, 3);
    OMEGA_H_CHECK(are_close(N_c, normalize(vector_3(1, 1, 0))));
  }
  {
    Few<Vector<3>, 3> N = {{1, 0, 0}, {1, 0, 0}, {1, 0, 0}};
    auto N_c = get_most_normal_normal(N, 3);
    OMEGA_H_CHECK(are_close(N_c, normalize(vector_3(1, 0, 0))));
  }
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
    OMEGA_H_CHECK(
        are_close(metric_desired_length(mi, q[i]), hi_expect[i], 1e-3));
  }
}

static void test_intersect_subset_metrics() {
  auto h1 = vector_2(1, 2);
  auto r1 = identity_matrix<2, 2>();
  auto h2 = vector_2(2, 3);
  auto r2 = rotate(PI / 4);
  auto m1 = compose_metric(r1, h1);
  auto m2 = compose_metric(r2, h2);
  OMEGA_H_CHECK(are_close(intersect_metrics(m2, m1), m1));
  OMEGA_H_CHECK(are_close(intersect_metrics(m1, m2), m1));
}

static void test_intersect_with_null() {
  auto q =
      rotate(PI / 4., vector_3(0, 0, 1)) * rotate(PI / 4., vector_3(0, 1, 0));
  auto m1 = compose_metric(q, vector_3(1, 1, 1e-3));
  auto m2 = zero_matrix<3, 3>();
  OMEGA_H_CHECK(are_close(intersect_metrics(m1, m2), m1));
  OMEGA_H_CHECK(are_close(intersect_metrics(m2, m1), m1));
}

static void test_intersect_degen_metrics() {
  test_intersect_with_null();
  // 2.a
  OMEGA_H_CHECK(are_close(intersect_metrics(diagonal(vector_3(1, 0, 0)),
                              diagonal(vector_3(0, 0, 1))),
      diagonal(vector_3(1, 0, 1))));
  // 2.b
  OMEGA_H_CHECK(are_close(intersect_metrics(diagonal(vector_3(1, 0, 0)),
                              diagonal(vector_3(2, 0, 0))),
      diagonal(vector_3(2, 0, 0))));
  // 3.a
  OMEGA_H_CHECK(are_close(intersect_metrics(diagonal(vector_3(1, 0, 0)),
                              diagonal(vector_3(2, 1, 0))),
      diagonal(vector_3(2, 1, 0))));
  // 3.b
  OMEGA_H_CHECK(are_close(intersect_metrics(diagonal(vector_3(1, 0, 0)),
                              diagonal(vector_3(0, 1, 2))),
      diagonal(vector_3(1, 1, 2))));
  // 4.a
  OMEGA_H_CHECK(are_close(intersect_metrics(diagonal(vector_3(1, 0, 2)),
                              diagonal(vector_3(2, 0, 1))),
      diagonal(vector_3(2, 0, 2))));
  // 4.b
  OMEGA_H_CHECK(are_close(intersect_metrics(diagonal(vector_3(1, 0, 2)),
                              diagonal(vector_3(2, 1, 0))),
      diagonal(vector_3(2, 1, 2))));
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
  test_intersect_degen_metrics();
}

static void test_interpolate_metrics() {
  auto a = repeat_symm(
      4, compose_metric(identity_matrix<2, 2>(), vector_2(1.0 / 100.0, 1.0)));
  auto b = repeat_symm(
      4, compose_metric(identity_matrix<2, 2>(), vector_2(1.0, 1.0)));
  auto c = interpolate_between_metrics(4, a, b, 0.0);
  OMEGA_H_CHECK(are_close(a, c));
  c = interpolate_between_metrics(4, a, b, 1.0);
  OMEGA_H_CHECK(are_close(b, c));
}

static void test_circumcenter() {
  Few<Vector<3>, 3> right_tri(
      {vector_3(0, 0, 0), vector_3(1, 0, 0), vector_3(0, 1, 0)});
  auto v0 = get_circumcenter_vector(simplex_basis<3, 2>(right_tri));
  OMEGA_H_CHECK(are_close(v0, vector_3(0.5, 0.5, 0)));
  Few<Vector<3>, 3> equal_tri(
      {vector_3(0, std::sqrt(3), 0), vector_3(-1, 0, 0), vector_3(1, 0, 0)});
  auto v1 = get_circumcenter_vector(simplex_basis<3, 2>(equal_tri));
  OMEGA_H_CHECK(are_close(v1, vector_3(0, -std::sqrt(3) * 2.0 / 3.0, 0)));
}

template <Int dim>
static void test_lie(Matrix<dim, dim> F) {
  auto log_F = log_polar(F);
  auto exp_log_F = exp_polar(log_F);
  OMEGA_H_CHECK(are_close(exp_log_F, F));
}

template <Int dim>
static Matrix<dim, dim> F_from_coords(Matrix<dim, dim + 1> new_simplex_coords) {
  Matrix<dim, dim + 1> old_simplex_basis_grads;
  for (Int i = 0; i < dim; ++i) {
    old_simplex_basis_grads(i, 0) = -1.0;
    for (Int j = 0; j < dim; ++j) {
      if (i == j)
        old_simplex_basis_grads(i, j + 1) = 1.0;
      else
        old_simplex_basis_grads(i, j + 1) = 0.0;
    }
  }
  return old_simplex_basis_grads * transpose(new_simplex_coords);
}

template <Int dim>
static void test_lie_F(Matrix<dim, dim + 1> new_simplex_coords) {
  test_lie(F_from_coords(new_simplex_coords));
}

static void test_lie() {
  test_lie(identity_matrix<1, 1>());
  test_lie(identity_matrix<2, 2>());
  test_lie(identity_matrix<3, 3>());
  test_lie_F<2>({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}});
  test_lie_F<2>({{0.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}});
  test_lie_F<2>({{0.0, 0.0}, {0.5, 0.0}, {0.0, 1.0}});
  test_lie_F<2>({{0.0, 0.0}, {1.0, 0.0}, {0.0, 2.0}});
  test_lie_F<2>({{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.5}});
  test_lie_F<2>({{0.0, 0.0}, {std::cos(1.0), std::sin(1.0)}, {0.0, 1.0}});
  test_lie_F<3>({{0.0, 0.0, 0.0}, {std::cos(1.0), std::sin(1.0), 0.0},
      {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
  test_lie_F<2>({{0.0, 0.0}, {std::cos(2.0), std::sin(2.0)}, {-1.0, 0.0}});
}

template <Int n>
static void test_positivize(Vector<n> pos) {
  auto neg = pos * -1.0;
  OMEGA_H_CHECK(are_close(positivize(pos), pos));
  OMEGA_H_CHECK(are_close(positivize(neg), pos));
}

static void test_positivize() {
  test_positivize(vector_3(1, 1, 1));
  test_positivize(vector_3(1, -1, 1));
  test_positivize(vector_2(-1, 1));
  test_positivize(vector_2(1, 1));
}

static void test_inball() {
  Few<Vector<1>, 2> regular_edge = {{-1.0}, {1.0}};
  auto inball1 = get_inball(regular_edge);
  OMEGA_H_CHECK(are_close(inball1.c, vector_1(0.0)));
  OMEGA_H_CHECK(are_close(inball1.r, 1.0));
  Few<Vector<2>, 3> regular_tri = {
      {-1.0, 0.0}, {1.0, 0.0}, {0.0, std::sqrt(3.0)}};
  auto inball2 = get_inball(regular_tri);
  OMEGA_H_CHECK(are_close(inball2.c, vector_2(0.0, std::sqrt(3.0) / 3.0)));
  OMEGA_H_CHECK(are_close(inball2.r, std::sqrt(3.0) / 3.0));
  Few<Vector<3>, 4> regular_tet = {{1, 0, -1.0 / std::sqrt(2.0)},
      {-1, 0, -1.0 / std::sqrt(2.0)}, {0, -1, 1.0 / std::sqrt(2.0)},
      {0, 1, 1.0 / std::sqrt(2.0)}};
  auto inball3 = get_inball(regular_tet);
  OMEGA_H_CHECK(are_close(inball3.c, vector_3(0.0, 0.0, 0.0)));
  OMEGA_H_CHECK(are_close(inball3.r, 2.0 / std::sqrt(24.0)));
}

static void test_volume_vert_gradients() {
  {
    Few<Vector<1>, 2> parent_edge = {{0.0}, {1.0}};
    auto grad0 = get_size_gradient(parent_edge, 0);
    OMEGA_H_CHECK(are_close(grad0, vector_1(-1.0)));
    auto grad1 = get_size_gradient(parent_edge, 1);
    OMEGA_H_CHECK(are_close(grad1, vector_1(1.0)));
  }
  {
    Few<Vector<2>, 3> parent_tri = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    auto grad0 = get_size_gradient(parent_tri, 0);
    OMEGA_H_CHECK(are_close(grad0, vector_2(-0.5, -0.5)));
    auto grad1 = get_size_gradient(parent_tri, 1);
    OMEGA_H_CHECK(are_close(grad1, vector_2(0.5, 0.0)));
    auto grad2 = get_size_gradient(parent_tri, 2);
    OMEGA_H_CHECK(are_close(grad2, vector_2(0.0, 0.5)));
  }
  {
    Few<Vector<3>, 4> parent_tet = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    auto os = 1.0 / 6.0;
    auto grad0 = get_size_gradient(parent_tet, 0);
    OMEGA_H_CHECK(are_close(grad0, vector_3(-os, -os, -os)));
    auto grad1 = get_size_gradient(parent_tet, 1);
    OMEGA_H_CHECK(are_close(grad1, vector_3(os, 0, 0)));
    auto grad2 = get_size_gradient(parent_tet, 2);
    OMEGA_H_CHECK(are_close(grad2, vector_3(0, os, 0)));
    auto grad3 = get_size_gradient(parent_tet, 3);
    OMEGA_H_CHECK(are_close(grad3, vector_3(0, 0, os)));
  }
}

template <Int dim>
static void test_svd_properties(Matrix<dim, dim> const A) {
  auto const svd = decompose_svd(A);
  OMEGA_H_CHECK(are_close(svd.U * svd.S * svd.V, A));
  OMEGA_H_CHECK(
      are_close(svd.U * transpose(svd.U), identity_matrix<dim, dim>()));
  OMEGA_H_CHECK(
      are_close(svd.V * transpose(svd.V), identity_matrix<dim, dim>()));
}

static void test_svd() {
  {
    auto const a = identity_matrix<1, 1>();
    auto const svd = decompose_svd(a);
    OMEGA_H_CHECK(are_close(svd.U(0, 0), 1.0));
    OMEGA_H_CHECK(are_close(svd.S(0, 0), 1.0));
    OMEGA_H_CHECK(are_close(svd.V(0, 0), 1.0));
  }
  {
    auto const a = identity_matrix<1, 1>() * 0.5;
    auto const svd = decompose_svd(a);
    OMEGA_H_CHECK(are_close(svd.U(0, 0), 1.0));
    OMEGA_H_CHECK(are_close(svd.S(0, 0), 0.5));
    OMEGA_H_CHECK(are_close(svd.V(0, 0), 1.0));
  }
  {
    auto const a = identity_matrix<2, 2>();
    auto const svd = decompose_svd(a);
    OMEGA_H_CHECK(are_close(svd.U, a));
    OMEGA_H_CHECK(are_close(svd.S, a));
    OMEGA_H_CHECK(are_close(svd.V, a));
  }
  {
    auto const I = identity_matrix<2, 2>();
    auto const a = I * 0.5;
    auto const svd = decompose_svd(a);
    OMEGA_H_CHECK(are_close(svd.U, I));
    OMEGA_H_CHECK(are_close(svd.S, a));
    OMEGA_H_CHECK(are_close(svd.V, I));
  }
  test_svd_properties(F_from_coords<2>({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}));
  test_svd_properties(F_from_coords<2>({{0.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}}));
  test_svd_properties(F_from_coords<2>({{0.0, 0.0}, {0.5, 0.0}, {0.0, 1.0}}));
  test_svd_properties(F_from_coords<2>({{0.0, 0.0}, {1.0, 0.0}, {0.0, 2.0}}));
  test_svd_properties(F_from_coords<2>({{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.5}}));
  test_svd_properties(F_from_coords<2>(
      {{0.0, 0.0}, {std::cos(1.0), std::sin(1.0)}, {0.0, 1.0}}));
  test_svd_properties(F_from_coords<3>({{0.0, 0.0, 0.0},
      {std::cos(1.0), std::sin(1.0), 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}));
  test_svd_properties(F_from_coords<2>(
      {{0.0, 0.0}, {std::cos(2.0), std::sin(2.0)}, {-1.0, 0.0}}));
}

static void test_quaternion(double angle, Vector<3> axis) {
  auto const tensor = rotate(angle, axis);
  auto const quaternion = quaternion_from_tensor(tensor);
  auto const tensor2 = tensor_from_quaternion(quaternion);
  auto const axis_angle = axis_angle_from_quaternion(quaternion);
  auto const angle2 = norm(axis_angle);
  auto axis2 = axis_angle;
  if (angle2 > DBL_EPSILON) axis2 /= angle2;
  OMEGA_H_CHECK(are_close(tensor, tensor2));
  OMEGA_H_CHECK(are_close(angle, angle2));
  if (angle > DBL_EPSILON) OMEGA_H_CHECK(are_close(axis, axis2));
  auto const quaternion2 = quaternion_from_axis_angle(angle * axis);
  OMEGA_H_CHECK(are_close(quaternion, quaternion2));
}

static void test_quaternions() {
  test_quaternion(0.0, vector_3(1.0, 0.0, 0.0));
  test_quaternion(PI, vector_3(1.0, 0.0, 0.0));
  test_quaternion(PI, vector_3(0.0, 1.0, 0.0));
  test_quaternion(PI, vector_3(0.0, 0.0, 1.0));
  test_quaternion(PI / 2.0, vector_3(1.0, 0.0, 0.0));
  test_quaternion(PI / 2.0, vector_3(0.0, 1.0, 0.0));
  test_quaternion(PI / 2.0, vector_3(0.0, 0.0, 1.0));
  test_quaternion(PI / 4.0, vector_3(1.0, 0.0, 0.0));
  test_quaternion(PI / 4.0, vector_3(0.0, 1.0, 0.0));
  test_quaternion(PI / 4.0, vector_3(0.0, 0.0, 1.0));
  test_quaternion(PI / 4.0, normalize(vector_3(1.0, 1.0, 1.0)));
  test_quaternion(PI / 4.0, normalize(vector_3(-1.0, 1.0, 1.0)));
  test_quaternion(PI / 4.0, normalize(vector_3(1.0, -1.0, 1.0)));
  test_quaternion(PI / 4.0, normalize(vector_3(-1.0, -1.0, 1.0)));
  test_quaternion(PI / 4.0, normalize(vector_3(1.0, 1.0, -1.0)));
  test_quaternion(PI / 4.0, normalize(vector_3(-1.0, 1.0, -1.0)));
  test_quaternion(PI / 4.0, normalize(vector_3(1.0, -1.0, -1.0)));
  test_quaternion(PI / 4.0, normalize(vector_3(-1.0, -1.0, -1.0)));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  test_edge_length();
  test_least_squares();
  test_power();
  test_cubic();
  test_form_ortho_basis();
  test_qr_decomps();
  test_eigen(matrix_1x1(42.0), vector_1(42.0));
  test_eigen_quadratic();
  test_eigen_cubic();
  test_eigen_jacobi();
  test_most_normal();
  test_intersect_metrics();
  test_interpolate_metrics();
  test_circumcenter();
  test_lie();
  test_positivize();
  test_inball();
  test_volume_vert_gradients();
  test_svd();
  test_quaternions();
}
