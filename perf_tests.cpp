#include <random>
#include "internal.hpp"

static UInt const nelems = 1;

static Reals random_reals(UInt n, Real from, Real to) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(from, to);
  HostReals ha(n);
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
  Reals alphas = random_reals(nelems, 0, PI / 2);
  Reals betas = random_reals(nelems, 0, PI / 2);
}

int main() {
  test_metric_decom();
}
