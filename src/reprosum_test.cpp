#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_array_ops.hpp"

static void test_repro_sum() {
  using namespace Omega_h;
  {
    Reals a = {0,1,2};
    Real res = repro_sum(a);
    printf("result %f\n", res);
    OMEGA_H_CHECK(are_close(res,3.0));
  }
  {
    Reals a({std::exp2(int(20)), std::exp2(int(-20))});
    Real sum = repro_sum(a);
    OMEGA_H_CHECK(sum == std::exp2(20) + std::exp2(int(-20)));
  }
  {
    const int n = 1'000'000;
    Write<Real> a(n);
    parallel_for(n, OMEGA_H_LAMBDA(int i) { a[i] = i; }, "setVals");
    auto const res = repro_sum(read(a));
    const double expected = (n-1)*static_cast<double>(n)/2.0;
    if(res != expected) fprintf(stderr, "expected %f != res %f\n", expected, res);
    OMEGA_H_CHECK(res == expected);
  }

}

int main(int argc, char** argv) {
  using namespace Omega_h;
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  test_repro_sum();
  return 0;
}
