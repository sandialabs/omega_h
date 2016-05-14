#include <random>
#include "internal.hpp"

static UInt const nelems = 100 * 1000;

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
    Matrix<3,3> l;
    decompose_eigen_qr(m, r, l);
    auto eigenv = max2(max2(l[0][0], l[1][1]), l[2][2]);
    write_eigenvs[i] = eigenv;
  };
  Now t0 = now();
  UInt niters = 30;
  for (UInt i = 0; i < niters; ++i)
    parallel_for(nelems, f1);
  Now t1 = now();
  std::cout << "eigendecomposition of " << nelems << " metric tensors "
    << niters << " times takes " << (t1 - t0) << " seconds\n";
  CHECK(are_close(Reals(write_eigenvs), Reals(nelems, square(anisotropy))));
}

static void test_repro_sum() {
  Reals inputs = random_reals(nelems, 0, 1e6);
  Real rsum;
  UInt niters = 500;
  Now t0 = now();
  for (UInt i = 0; i < niters; ++i)
    rsum = repro_sum(inputs);
  Now t1 = now();
  std::cout << "adding " << nelems << " reals " << niters << " times "
    << "takes " << (t1 - t0) << " seconds\n";
}

int main() {
  test_metric_decom();
  test_repro_sum();
}
