#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_sort.hpp"

int main(int argc, char** argv) {
  using namespace Omega_h;
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  {
    Reals a = {1,0,2};
    bool res = is_sorted(a);
    OMEGA_H_CHECK(res == false);
  }
  {
    Reals a = {1,0,1,2};
    LO res = find_last(a,1.0);
    OMEGA_H_CHECK(res == 2);
  }
  {
    Reals a = {1,0.0,1,2};
    Reals b = {1,0.1,1,2};
    bool res = are_close_abs(a,b,1e-3);
    OMEGA_H_CHECK(res == false);
    Reals c = {1,0.000001,1,2};
    res = are_close_abs(a,c,1e-3);
    OMEGA_H_CHECK(res == true);
  }
  {
    const auto tol = 1e-3;
    const auto floor = 1e-9;
    Reals a = {1,0.0,1,2};
    Reals b = {1,0.1,1,2};
    bool res = are_close(a,b,tol,floor);
    OMEGA_H_CHECK(res == false);
    const auto small = 0.000001;
    Reals c = {1,small,1,2};
    res = are_close(a,c,tol,floor);
    const auto expected = are_close(0.1,small,tol,floor);
    OMEGA_H_CHECK(res == expected);
  }
  {
    Reals a = {1,0.1,1,2};
    auto const res = get_max(a);
    OMEGA_H_CHECK(res == 2);
  }
  {
    Reals a = {1,0.1,1,2};
    auto const res = get_min(a);
    OMEGA_H_CHECK(res == 0.1);
  }
  {
    const int n = 1000;
    Write<Real> a(n);
    parallel_for(n, OMEGA_H_LAMBDA(int i) { a[i] = i; }, "setVals");
    a.set(42,0.1);
    auto const res = get_min(read(a));
    OMEGA_H_CHECK(res == 0.1);
  }
  {
    Reals a = {1,0.1,1,2};
    auto const res = get_sum(a);
    OMEGA_H_CHECK(res == 4.1);
  }
  {
    const int n = 1000;
    Write<LO> a(n);
    parallel_for(n, OMEGA_H_LAMBDA(int i) { a[i] = i; }, "setVals");
    auto const res = get_sum(read(a));
    const int expected = (n-1)*(n)/2;
    if( expected != res )
      fprintf(stderr, "expected %d res %d\n", expected, res);
    OMEGA_H_CHECK(res == expected);
  }
  {
    LOs a = {1,5,1,2};
    LOs b = {1,5,1,2};
    LOs c = {1,5,2,2};
    auto res = (a==b);
    OMEGA_H_CHECK(res == true);
    res = (a==c);
    OMEGA_H_CHECK(res == false);
  }
  {
    Write<LO> a = {1,2,0};
    LOs b = {1,2,3};
    a.set(2,3);
    auto res = (read(a)==b);
    OMEGA_H_CHECK(res == true);
    auto val = a.get(1);
    OMEGA_H_CHECK(val == 2);
  }
  {
    Read<I32> a = {45,2,0,3};
    LOs expectedPerm = {3,1,0,2};
    LOs expectedFan = {0,1,2,3,4};
    LOs expectedUniq = {0,2,3,45};
    LOs perm, fan;
    Read<I32> uniq;
    sort_small_range(a,&perm,&fan,&uniq);
    OMEGA_H_CHECK(perm == expectedPerm);
    OMEGA_H_CHECK(fan == expectedFan);
    OMEGA_H_CHECK(uniq == expectedUniq);
  }
  return 0;
}
