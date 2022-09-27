#include "Omega_h_adj.hpp"
#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_expr.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_sort.hpp"
#include "Omega_h_atomics.hpp"

using namespace Omega_h;

static void test_write() {
  auto w = Write<Real>(100, "foo");
  OMEGA_H_CHECK(w.size() == 100);
  OMEGA_H_CHECK(w.name() == "foo");
}

static void test_int128() {
  Int128 a(INT64_MAX);
  auto b = a + a;
  b = b + b;
  b = b + b;
  b = b >> 3;
  OMEGA_H_CHECK(b == a);
}

static void test_atomic() {
  const auto a2b = LOs({0,1,1,3,2});
  const auto na = a2b.size();
  Write<LO> degrees(4, 0);
  auto count = OMEGA_H_LAMBDA(LO a) {
    auto const b = a2b[a];
    atomic_increment(&degrees[b]);
  };
  parallel_for(na, std::move(count));
  OMEGA_H_CHECK(read(degrees) == LOs({1,2,1,1}));
}

static void test_repro_sum() {
  Reals a({std::exp2(int(20)), std::exp2(int(-20))});
  Real sum = repro_sum(a);
  OMEGA_H_CHECK(sum == std::exp2(20) + std::exp2(int(-20)));
}

static void test_sort() {
  {
    LOs a({0, 1});
    LOs perm = sort_by_keys(a);
    OMEGA_H_CHECK(perm == LOs({0, 1}));
  }
  {
    LOs a({0, 2, 0, 1});
    LOs perm = sort_by_keys(a, 2);
    OMEGA_H_CHECK(perm == LOs({1, 0}));
  }
  {
    LOs a({0, 2, 1, 1});
    LOs perm = sort_by_keys(a, 2);
    OMEGA_H_CHECK(perm == LOs({0, 1}));
  }
  {
    LOs a({1, 2, 3, 1, 2, 2, 3, 0, 0});
    LOs perm = sort_by_keys(a, 3);
    OMEGA_H_CHECK(perm == LOs({1, 0, 2}));
  }
}

static void test_sort_small_range() {
  Read<I32> in({10, 100, 1000, 10, 100, 1000, 10, 100, 1000});
  LOs perm;
  LOs fan;
  Read<I32> uniq;
  fprintf(stderr, "test_sort_small 0.1\n");
  sort_small_range(in, &perm, &fan, &uniq);
  fprintf(stderr, "test_sort_small 0.2\n");
  OMEGA_H_CHECK(perm == LOs({0, 3, 6, 1, 4, 7, 2, 5, 8}));
  OMEGA_H_CHECK(fan == LOs({0, 3, 6, 9}));
  OMEGA_H_CHECK(uniq == Read<I32>({10, 100, 1000}));
  in = Read<I32>({});
  sort_small_range(in, &perm, &fan, &uniq);
  OMEGA_H_CHECK(perm == LOs({}));
  OMEGA_H_CHECK(fan == LOs({0}));
  OMEGA_H_CHECK(uniq == Read<I32>({}));
}

static void test_scan() {
  {
    LOs scanned = offset_scan(LOs(3, 1));
    OMEGA_H_CHECK(scanned == Read<LO>(4, 0, 1));
  }
  {
    LOs scanned = offset_scan(Read<I8>(3, 1));
    OMEGA_H_CHECK(scanned == Read<LO>(4, 0, 1));
  }
}

static void test_fan_and_funnel() {
  OMEGA_H_CHECK(invert_funnel(LOs({0, 0, 1, 1, 2, 2}), 3) == LOs({0, 2, 4, 6}));
  OMEGA_H_CHECK(invert_fan(LOs({0, 2, 4, 6})) == LOs({0, 0, 1, 1, 2, 2}));
  OMEGA_H_CHECK(invert_funnel(LOs({0, 0, 0, 2, 2, 2}), 3) == LOs({0, 3, 3, 6}));
  OMEGA_H_CHECK(invert_fan(LOs({0, 3, 3, 6})) == LOs({0, 0, 0, 2, 2, 2}));
  OMEGA_H_CHECK(invert_funnel(LOs({0, 0, 0, 0, 0, 0}), 3) == LOs({0, 6, 6, 6}));
  OMEGA_H_CHECK(invert_fan(LOs({0, 6, 6, 6})) == LOs({0, 0, 0, 0, 0, 0}));
  OMEGA_H_CHECK(invert_funnel(LOs({2, 2, 2, 2, 2, 2}), 3) == LOs({0, 0, 0, 6}));
  OMEGA_H_CHECK(invert_fan(LOs({0, 0, 0, 6})) == LOs({2, 2, 2, 2, 2, 2}));
}

static void test_permute() {
  Reals data({0.1, 0.2, 0.3, 0.4});
  LOs perm({3, 2, 1, 0});
  Reals permuted = unmap(perm, data, 1);
  OMEGA_H_CHECK(permuted == Reals({0.4, 0.3, 0.2, 0.1}));
  Reals back = permute(permuted, perm, 1);
  OMEGA_H_CHECK(back == data);
}

static void test_invert_map() {
  {
    LOs hl2l({}, "hl2l");
    auto l2hl = invert_map_by_atomics(hl2l, 4);
    OMEGA_H_CHECK(l2hl.a2ab == LOs(5, 0));
    OMEGA_H_CHECK(l2hl.ab2b == LOs({}));
  }
  {
    LOs hl2l({0, 1, 2, 3}, "hl2l");
    auto l2hl = invert_map_by_atomics(hl2l, 4);
    OMEGA_H_CHECK(l2hl.a2ab == LOs(5, 0, 1));
    OMEGA_H_CHECK(l2hl.ab2b == LOs(4, 0, 1));
  }
  {
    auto l2hl = invert_map_by_atomics(LOs({0,1,2,2,3,0}),4);
    OMEGA_H_CHECK(l2hl.a2ab == LOs({0, 2, 3, 5, 6}));
    OMEGA_H_CHECK(l2hl.ab2b == LOs({0, 5, 1, 2, 3, 4}));
  }
  {
    LOs hl2l({}, "hl2l");
    auto l2hl = invert_map_by_sorting(hl2l, 4);
    OMEGA_H_CHECK(l2hl.a2ab == LOs(5, 0));
    OMEGA_H_CHECK(l2hl.ab2b == LOs({}));
  }
  {
    LOs hl2l({0, 1, 2, 3}, "hl2l");
    auto l2hl = invert_map_by_sorting(hl2l, 4);
    OMEGA_H_CHECK(l2hl.a2ab == LOs(5, 0, 1));
    OMEGA_H_CHECK(l2hl.ab2b == LOs(4, 0, 1));
  }
  {
    LOs hl2l({1, 0, 1, 0}, "hl2l");
    auto l2hl = invert_map_by_sorting(hl2l, 2);
    OMEGA_H_CHECK(l2hl.a2ab == LOs({0, 2, 4}));
    OMEGA_H_CHECK(l2hl.ab2b == LOs({1, 3, 0, 2}));
  }
}

static void test_invert_adj() {
  Adj tris2verts(LOs({0, 1, 2, 2, 3, 0}));
  Adj verts2tris = invert_adj(tris2verts, 3, 4, 2, 0);
  OMEGA_H_CHECK(verts2tris.a2ab == offset_scan(LOs({2, 1, 2, 1})));
  OMEGA_H_CHECK(verts2tris.ab2b == LOs({0, 1, 0, 0, 1, 1}));
  OMEGA_H_CHECK(
      verts2tris.codes ==
      Read<I8>({make_code(0, 0, 0), make_code(0, 0, 2), make_code(0, 0, 1),
          make_code(0, 0, 2), make_code(0, 0, 0), make_code(0, 0, 1)}));
}

static void test_injective_map() {
  LOs primes2ints({2, 3, 5, 7});
  LOs ints2primes = invert_injective_map(primes2ints, 8);
  OMEGA_H_CHECK(ints2primes == LOs({-1, -1, 0, 1, -1, 2, -1, 3}));
}

void test_binary_search(LOs a, LO val, LO expect) {
  auto size = a.size();
  auto f = OMEGA_H_LAMBDA(LO) {
    auto res = binary_search(a, val, size);
    OMEGA_H_CHECK(res == expect);
  };
  parallel_for(1, f);
}

static void test_binary_search() {
  auto a = LOs({20, 30, 40, 50, 90, 100});
  for (LO i = 0; i < a.size(); ++i) test_binary_search(a, a.get(i), i);
  test_binary_search(a, 44, -1);
  test_binary_search(a, 15, -1);
  test_binary_search(a, 10000, -1);
}

static void test_is_sorted() {
  OMEGA_H_CHECK(is_sorted(LOs({})));
  OMEGA_H_CHECK(is_sorted(Reals({42.0})));
  OMEGA_H_CHECK(is_sorted(LOs({0, 1, 2})));
  OMEGA_H_CHECK(!is_sorted(Reals({0.2, 0.1, 0.3, 0.4})));
  OMEGA_H_CHECK(is_sorted(Reals({0.1, 0.1, 0.1, 0.1})));
}

static void test_linpart() {
  GO total = 7;
  I32 comm_size = 2;
  OMEGA_H_CHECK(linear_partition_size(total, comm_size, 0) == 4);
  OMEGA_H_CHECK(linear_partition_size(total, comm_size, 1) == 3);
  Read<GO> globals({6, 5, 4, 3, 2, 1, 0});
  auto remotes = globals_to_linear_owners(globals, total, comm_size);
  OMEGA_H_CHECK(remotes.ranks == Read<I32>({1, 1, 1, 0, 0, 0, 0}));
  OMEGA_H_CHECK(remotes.idxs == Read<I32>({2, 1, 0, 3, 2, 1, 0}));
}

static void test_expand() {
  auto fan = offset_scan(LOs({2, 1, 3}));
  Reals data({2.2, 3.14, 42.0});
  OMEGA_H_CHECK(
      expand(data, fan, 1) == Reals({2.2, 2.2, 3.14, 42.0, 42.0, 42.0}));
}

static void test_find_last() {
  auto a = LOs({0, 3, 55, 12});
  OMEGA_H_CHECK(find_last(a, 98) < 0);
  OMEGA_H_CHECK(find_last(a, 12) == 3);
  OMEGA_H_CHECK(find_last(a, 55) == 2);
  OMEGA_H_CHECK(find_last(a, 3) == 1);
  OMEGA_H_CHECK(find_last(a, 0) == 0);
}

static void test_scalar_ptr() {
  Vector<2> v;
  OMEGA_H_CHECK(scalar_ptr(v) == &v[0]);
  auto const& v2 = v;
  OMEGA_H_CHECK(scalar_ptr(v2) == &v2[0]);
  OMEGA_H_CHECK(scalar_ptr(v2) == &v2(0));
  Matrix<5, 4> m;
  OMEGA_H_CHECK(scalar_ptr(m) == &m[0][0]);
  auto const& m2 = m;
  OMEGA_H_CHECK(scalar_ptr(m2) == &m2[0][0]);
  OMEGA_H_CHECK(scalar_ptr(m2) == &m2(0, 0));
}

static void test_expr() {
  using Omega_h::any;
  using Omega_h::any_cast;
  ExprReader reader(4, 3);
  OMEGA_H_CHECK(any_cast<Real>(reader.read_string("1.0", "test0")) == 1.0);
  OMEGA_H_CHECK(any_cast<Real>(reader.read_string("1 + 1", "test1")) == 2.0);
  reader.register_variable("pi", any(Real(3.14159)));
  OMEGA_H_CHECK(any_cast<Real>(reader.read_string("pi", "test2")) == 3.14159);
  reader.register_variable("j", any(vector_3(0, 1, 0)));
  OMEGA_H_CHECK(
      are_close(any_cast<Vector<3>>(reader.read_string("pi * j", "test3")),
          vector_3(0, 3.14159, 0)));
  reader.register_variable("x", any(Reals({0, 1, 2, 3})));
  OMEGA_H_CHECK(are_close(any_cast<Reals>(reader.read_string("x^2", "test4")),
      Reals({0, 1, 4, 9})));
  reader.register_variable(
      "v", any(Reals({0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0})));
  OMEGA_H_CHECK(are_close(any_cast<Reals>(reader.read_string("v(1)", "test5")),
      Reals({0, 1, 2, 3})));
  OMEGA_H_CHECK(
      are_close(any_cast<Reals>(reader.read_string("v - 1.5 * j", "test6")),
          Reals({0, -1.5, 0, 0, -0.5, 0, 0, 0.5, 0, 0, 1.5, 0})));
  OMEGA_H_CHECK(are_close(
      any_cast<Vector<3>>(reader.read_string("vector(0, 1, 2)", "test7")),
      vector_3(0, 1, 2)));
  OMEGA_H_CHECK(
      are_close(any_cast<Reals>(reader.read_string("vector(x, 0, 0)", "test8")),
          Reals({0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0})));
  OMEGA_H_CHECK(
      are_close(any_cast<Real>(reader.read_string("exp(0)", "test9")), 1.0));
  OMEGA_H_CHECK(
      are_close(any_cast<Reals>(reader.read_string("exp(x)", "test10")),
          Reals({1.0, std::exp(1.0), std::exp(2.0), std::exp(3.0)})));
}

static Omega_h::any test_expr2(
    ExprEnv& env, std::string const& expr, std::string const& test_name) {
  ExprOpsReader reader;
  auto op = any_cast<OpPtr>(reader.read_string(expr, test_name));
  return op->eval(env);
}

static void test_expr2() {
  using Omega_h::any;
  using Omega_h::any_cast;
  ExprEnv env(4, 3);
  OMEGA_H_CHECK(any_cast<Real>(test_expr2(env, "1.0", "test0")) == 1.0);
  OMEGA_H_CHECK(any_cast<Real>(test_expr2(env, "1 + 1", "test1")) == 2.0);
  env.register_variable("pi", any(Real(3.14159)));
  OMEGA_H_CHECK(any_cast<Real>(test_expr2(env, "pi", "test2")) == 3.14159);
  env.register_variable("j", any(vector_3(0, 1, 0)));
  OMEGA_H_CHECK(
      are_close(any_cast<Vector<3>>(test_expr2(env, "pi * j", "test3")),
          vector_3(0, 3.14159, 0)));
  env.register_variable("x", any(Reals({0, 1, 2, 3})));
  OMEGA_H_CHECK(are_close(
      any_cast<Reals>(test_expr2(env, "x^2", "test4")), Reals({0, 1, 4, 9})));
  env.register_variable("v", any(Reals({0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0})));
  OMEGA_H_CHECK(are_close(
      any_cast<Reals>(test_expr2(env, "v(1)", "test5")), Reals({0, 1, 2, 3})));
  OMEGA_H_CHECK(
      are_close(any_cast<Reals>(test_expr2(env, "v - 1.5 * j", "test6")),
          Reals({0, -1.5, 0, 0, -0.5, 0, 0, 0.5, 0, 0, 1.5, 0})));
  OMEGA_H_CHECK(are_close(
      any_cast<Vector<3>>(test_expr2(env, "vector(0, 1, 2)", "test7")),
      vector_3(0, 1, 2)));
  OMEGA_H_CHECK(
      are_close(any_cast<Reals>(test_expr2(env, "vector(x, 0, 0)", "test8")),
          Reals({0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0})));
  OMEGA_H_CHECK(
      are_close(any_cast<Real>(test_expr2(env, "exp(0)", "test9")), 1.0));
  OMEGA_H_CHECK(are_close(any_cast<Reals>(test_expr2(env, "exp(x)", "test10")),
      Reals({1.0, std::exp(1.0), std::exp(2.0), std::exp(3.0)})));
}

static void test_array_from_kokkos() {
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::View<double**> managed(
      Kokkos::ViewAllocateWithoutInitializing("view"), 10, 10);
  Kokkos::View<double*> unmanaged(managed.data(), managed.span());
  Omega_h::Write<double> unmanaged_array(unmanaged);
  OMEGA_H_CHECK(unmanaged_array.exists());
  Kokkos::View<double*> zero_span("zero_span", 0);
  Omega_h::Write<double> zero_span_array(zero_span);
  OMEGA_H_CHECK(zero_span_array.exists());
  Kokkos::View<double*> uninitialized;
  Omega_h::Write<double> uninitialized_array(uninitialized);
  OMEGA_H_CHECK(!uninitialized_array.exists());
#endif
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  test_write();
  test_atomic();
  test_int128();
  test_repro_sum();
  test_sort();
  fprintf(stderr, "0.1\n");
  test_sort_small_range();
  fprintf(stderr, "0.2\n");
  test_scan();
  fprintf(stderr, "0.3\n");
  test_fan_and_funnel();
  test_permute();
  test_invert_map();
  test_invert_adj();
  test_injective_map();
  test_binary_search();
  test_is_sorted();
  test_linpart();
  test_expand();
  test_find_last();
  test_scalar_ptr();
  test_expr();
  test_expr2();
  test_array_from_kokkos();
}
