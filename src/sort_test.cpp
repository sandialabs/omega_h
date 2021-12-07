#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include "Omega_h_library.hpp"
#include "Omega_h_array_ops.hpp"
//#include "Omega_h_sort.hpp"

int main(int argc, char** argv) {
  using namespace Omega_h;
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  {
    //Write<LO> perm(n, 0, 1);

    auto space = Kokkos::Experimental::SYCL();
    const auto q = *space.impl_internal_space_instance()->m_queue;

    using namespace oneapi::dpl::execution;
    using namespace sycl;
    auto policy = make_device_policy(q);

    Write<LO> perm = {43, 0, 3};
    LOs expected = {0, 3, 43};
    oneapi::dpl::sort(
        policy, perm.data(), perm.data() + perm.size(),
        [](auto lhs, auto rhs) { return lhs < rhs; });

    HostRead<LO> h_perm(read(perm));
    for(int i=0; i<h_perm.size(); i++)
      printf("%d ", h_perm[i]);
    printf("\n");

    Read<LO> r_expected(expected);
    Read<LO> r_perm(perm);
    OMEGA_H_CHECK(r_perm == r_expected);
  }
//  }
//  {
//    LOs a({0, 2, 0, 1});
//    LOs perm = sort_by_keys(a, 2);
//    OMEGA_H_CHECK(perm == LOs({1, 0}));
//  }
//  {
//    LOs a({0, 2, 1, 1});
//    LOs perm = sort_by_keys(a, 2);
//    OMEGA_H_CHECK(perm == LOs({0, 1}));
//  }
//  {
//    LOs a({1, 2, 3, 1, 2, 2, 3, 0, 0});
//    LOs perm = sort_by_keys(a, 3);
//    OMEGA_H_CHECK(perm == LOs({1, 0, 2}));
//  }
  return 0;
}
