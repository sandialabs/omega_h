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
    Reals a = {1,0.1,1,2};
    auto const res = get_sum(a);
    OMEGA_H_CHECK(res == 4.1);
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
    Read<I32> a = {1,2,0,3};
    LOs perm, fan;
    Read<I32> b;
    //sort_small_range(a,&perm,&fan,&b);
    //OMEGA_H_CHECK(res == 4.1);
  }
  { //the following works, but the call to number_same_values in sort.cpp fails
    Read<I32> a = {1,1,0,1};
    Write<LO> perm(a.size()+1,0);
    LOs expected = {0,1,2,2,3};
    I32 value = 1;
    auto transform = OMEGA_H_LAMBDA(LO i)->LO {
      return a[i] == value ? LO(1) : LO(0);
    };
    Kokkos::parallel_scan(
      Kokkos::RangePolicy<>(0, a.size() ),
      KOKKOS_LAMBDA(int i, LO& update, const bool final) {
        update += transform(i);
        if(final) perm[i+1] = update;
      });
    HostRead<LO> h_perm(perm);
    for(int i=0; i<h_perm.size(); i++)
      printf("%d ", h_perm[i]);
    printf("\n");
    OMEGA_H_CHECK(read(perm) == expected);
  }
  return 0;
}
