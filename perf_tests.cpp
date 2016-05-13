#include <random>
#include "internal.hpp"

static UInt const nelems = 1;

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
  Write<Real> write_eigenv(nelems);
  auto f1 = LAMBDA(UInt i) {
    auto m = get_symm<3>(metrics, i);
    Matrix<3,3> r;
    Matrix<3,3> l;
    decompose_eigen_qr(m, r, l);
    write_eigenv[i] = max2(max2(l[0][0], l[1][1]), l[2][2]);
  };
  parallel_for(nelems, f1);
}

int main() {
  test_metric_decom();
}
