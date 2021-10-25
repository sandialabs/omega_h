#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_int128.hpp"


namespace sample {
struct Int128Wrap {
  Omega_h::Int128 i128;

  KOKKOS_INLINE_FUNCTION   // Default constructor - Initialize to 0
  Int128Wrap() {
    i128 = Omega_h::Int128(0);
  }
  KOKKOS_INLINE_FUNCTION   // Copy Constructor
  Int128Wrap(const Int128Wrap & rhs) {
    i128 = rhs.i128;
  }
  KOKKOS_INLINE_FUNCTION   // add operator
  Int128Wrap& operator += (const Int128Wrap& src) {
    i128 = i128 + src.i128;
    return *this;
  }
  KOKKOS_INLINE_FUNCTION   // volatile add operator
  void operator += (const volatile Int128Wrap& src) volatile {
    i128 = i128 + src.i128;
  }
};
}

namespace Kokkos { //reduction identity must be defined in Kokkos namespace
template<>
struct reduction_identity< sample::Int128Wrap > {
   KOKKOS_FORCEINLINE_FUNCTION static sample::Int128Wrap sum() {
      return sample::Int128Wrap();
   }
};
}

int main(int argc, char** argv) {
  using namespace Omega_h;
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Reals a = {0,1,2};

  sample::Int128Wrap res;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<>(0, a.size() ),
    KOKKOS_LAMBDA(int i, sample::Int128Wrap& update) {
      update.i128 = Int128::from_double(a[i],1.0);
    }, Kokkos::Sum< sample::Int128Wrap >(res) );
   printf("result %zu %zu\n", res.i128.high, res.i128.low);

  return 0;
}
