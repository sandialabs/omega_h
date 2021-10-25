#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_bbox.hpp"
#include "Omega_h_bbox.hpp"

using namespace Omega_h;

namespace sample {
struct Int128Wrap {
  Int128 i128;

  KOKKOS_INLINE_FUNCTION   // Default constructor - Initialize to 0
  Int128Wrap() {
    i128 = Int128(0);
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
    const auto foo = i128 + src.i128; //FIXME does not compile, complains about no matching ctor for 'i128'
    i128 = foo;
  }
};
}

namespace Kokkos { //reduction identity must be defined in Kokkos namespace
template<>
struct reduction_identity< Omega_h::Int128 > {
   KOKKOS_FORCEINLINE_FUNCTION static Omega_h::Int128 sum() {
      return Omega_h::Int128();
   }
};
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Int128 a = {0,1,3};

  sample::Int128Wrap res;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<>(0, a.size() ),
    KOKKOS_LAMBDA(int i, sample::Int128Wrap& update) {
      update.i128 = transform(i);
    }, Kokkos::Sum< sample::Int128Wrap >(res) );
  return res.i128;

  return 0;
}
