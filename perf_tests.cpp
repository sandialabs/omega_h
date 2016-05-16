#include <random>
#include <cstdio>
#include "internal.hpp"

static Int const nelems = 1000 * 1000;

static Reals random_reals(Int n, Real from, Real to) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(from, to);
  HostWrite<Real> ha(n);
  for (Int i = 0; i < n; ++i)
    ha[i] = dis(gen);
  return ha.write();
}

static void test_metric_decom() {
/* we'll generate metric tensors which have random
   rotation matrices (choose two random Euler angles.
   this is not uniform over the sphere, but thats not important),
   and desired lengths (1,1,a^2), where (a) is the desired
   amount of anisotropy (a=1000 implies 1000:1 ratio). */
  double anisotropy = 1000;
  Reals alphas = random_reals(nelems, 0, PI / 2);
  Reals betas = random_reals(nelems, 0, PI / 2);
  Write<Real> write_metrics(nelems * 6);
  auto f0 = LAMBDA(Int i) {
    auto r = rotate(alphas[i], vector_3(0, 0, 1)) *
             rotate( betas[i], vector_3(0, 1, 0));
    auto l = matrix_3x3(
        1, 0, 0,
        0, 1, 0,
        0, 0, square(anisotropy));
    auto m = r * l * transpose(r);
    set_symm(write_metrics, i, m);
  };
  parallel_for(nelems, f0);
  Reals metrics(write_metrics);
  /* now, decompose the metrics and get the largest
     eigenvalue of each */
  Write<Real> write_eigenvs(nelems);
  auto f1 = LAMBDA(Int i) {
    auto m = get_symm<3>(metrics, i);
    Matrix<3,3> r;
    Vector<3> l;
    bool ok = decompose_eigen(m, r, l);
    CHECK(ok);
    auto eigenv = max2(max2(l[0], l[1]), l[2]);
    write_eigenvs[i] = eigenv;
  };
  Now t0 = now();
  Int niters = 3;
  for (Int i = 0; i < niters; ++i)
    parallel_for(nelems, f1);
  Now t1 = now();
  std::cout << "eigendecomposition of " << nelems << " metric tensors "
    << niters << " times takes " << (t1 - t0) << " seconds\n";
  CHECK(are_close(Reals(write_eigenvs), Reals(nelems, square(anisotropy))));
}

unsigned uniform(unsigned m); /* Returns a random integer 0 <= uniform(m) <= m-1 */

/* Fisher-Yates shuffle, a.k.a Knuth shuffle */
static Read<Int> random_perm(Int n)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> dis;
  /* start with identity permutation */
  HostWrite<Int> permutation(make_linear<Int>(n, 0, 1));
  for (Int i = 0; i + 1 < n; i++) {
    Int j = dis(gen) % (n - i); /* A random integer such that 0 ≤ j < n-i*/
    std::swap(permutation[i], permutation[i + j]);
  }
  return Read<Int>(permutation.write());
}

static void test_repro_sum() {
  Reals inputs = random_reals(nelems, 0, 1e100);
  Real rs = 0, s = 0;
  Int niters = 100;
  {
    Now t0 = now();
    for (Int i = 0; i < niters; ++i)
      rs = repro_sum(inputs);
    Now t1 = now();
    std::cout << "reproducibly adding " << nelems << " reals " << niters << " times "
      << "takes " << (t1 - t0) << " seconds\n";
  }
  {
    Now t0 = now();
    for (Int i = 0; i < niters; ++i)
      s = sum(inputs);
    Now t1 = now();
    std::cout << "adding " << nelems << " reals " << niters << " times "
      << "takes " << (t1 - t0) << " seconds\n";
  }
  CHECK(are_close(s, rs));
  Read<Int> p = random_perm(nelems);
  Write<Real> write_shuffled(nelems);
  auto f = LAMBDA(Int i) {
    write_shuffled[i] = inputs[p[i]];
  };
  parallel_for(nelems, f);
  Reals shuffled(write_shuffled);
  Real rs2 = repro_sum(shuffled);
  Real s2 = sum(shuffled);
  CHECK(are_close(s2, rs2));
  CHECK(rs == rs2); /* bitwise reproducibility ! */
  if (s == s2)
    std::cerr << "warning: the naive sum gave the same answer\n";
}

template <typename T>
static Read<T> random_ints(Int n, T from, T to) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<T> dis(from, to);
  HostWrite<T> ha(n);
  for (Int i = 0; i < n; ++i)
    ha[i] = dis(gen);
  return ha.write();
}

template <Int N>
static void test_sort_n() {
  LOs a = random_ints<LO>(nelems * N, 0, nelems);
  LOs perm;
  Int niters = 5;
  Now t0 = now();
  for (Int i = 0; i < niters; ++i)
    perm = sort_by_keys<LO,N>(a);
  Now t1 = now();
  std::cout << "sorting " << nelems << " sets of "
    << N << " integers " << niters
    << " times takes " << (t1 - t0) << " seconds\n";
}

static void test_sort() {
  test_sort_n<1>();
  test_sort_n<2>();
  test_sort_n<3>();
}

int main(int argc, char** argv) {
  init(argc, argv);
  test_metric_decom();
  test_repro_sum();
  test_sort();
  fini();
}
