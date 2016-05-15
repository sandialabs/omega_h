#include <random>
#include <cstdio>
#include "internal.hpp"

static UInt const nelems = 1000 * 1000;

static Reals random_reals(UInt n, Real from, Real to) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(from, to);
  HostWrite<Real> ha(n);
  for (UInt i = 0; i < n; ++i)
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
  auto f0 = LAMBDA(UInt i) {
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
  auto f1 = LAMBDA(UInt i) {
    auto m = get_symm<3>(metrics, i);
    Matrix<3,3> r;
    Vector<3> l;
    decompose_eigen_cubic(m, r, l);
    auto eigenv = max2(max2(l[0], l[1]), l[2]);
    write_eigenvs[i] = eigenv;
  };
  Now t0 = now();
  UInt niters = 3;
  for (UInt i = 0; i < niters; ++i)
    parallel_for(nelems, f1);
  Now t1 = now();
  std::cout << "eigendecomposition of " << nelems << " metric tensors "
    << niters << " times takes " << (t1 - t0) << " seconds\n";
  CHECK(are_close(Reals(write_eigenvs), Reals(nelems, square(anisotropy))));
}

unsigned uniform(unsigned m); /* Returns a random integer 0 <= uniform(m) <= m-1 */

/* Fisher-Yates shuffle, a.k.a Knuth shuffle */
static Read<UInt> random_perm(UInt n)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<UInt> dis;
  /* start with identity permutation */
  HostWrite<UInt> permutation(make_linear<UInt>(n, 0, 1));
  for (UInt i = 0; i <= n-2; i++) {
    UInt j = dis(gen) % (n - i); /* A random integer such that 0 ≤ j < n-i*/
    std::swap(permutation[i], permutation[i + j]);
  }
  return Read<UInt>(permutation.write());
}

static void test_repro_sum() {
  Reals inputs = random_reals(nelems, 0, 1e6);
  Real rs = 0, s = 0;
  UInt niters = 50;
  {
    Now t0 = now();
    for (UInt i = 0; i < niters; ++i)
      rs = repro_sum(inputs);
    Now t1 = now();
    std::cout << "reproducibly adding " << nelems << " reals " << niters << " times "
      << "takes " << (t1 - t0) << " seconds\n";
  }
  {
    Now t0 = now();
    for (UInt i = 0; i < niters; ++i)
      s = sum(inputs);
    Now t1 = now();
    std::cout << "adding " << nelems << " reals " << niters << " times "
      << "takes " << (t1 - t0) << " seconds\n";
  }
  CHECK(are_close(s, rs));
  Read<UInt> p = random_perm(nelems);
  Write<Real> write_shuffled(nelems);
  auto f = LAMBDA(UInt i) {
    write_shuffled[i] = inputs[p[i]];
  };
  parallel_for(nelems, f);
  Reals shuffled(write_shuffled);
  Real rs2 = repro_sum(shuffled);
  Real s2 = sum(shuffled);
  CHECK(are_close(s2, rs2));
  CHECK(s != s2);
  CHECK(rs == rs2); /* bitwise reproducibility ! */
}

int main(int argc, char** argv) {
  init(argc, argv);
  test_metric_decom();
  test_repro_sum();
  fini();
}
