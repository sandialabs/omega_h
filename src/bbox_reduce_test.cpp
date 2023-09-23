#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_bbox.hpp"
#include "Omega_h_bbox.hpp"

using namespace Omega_h;

template <Int dim>
struct GetBBoxOp2 {
  Reals coords;
  GetBBoxOp2(Reals coords_in) : coords(coords_in) {}
  OMEGA_H_INLINE BBox<dim> operator()(LO i) const {
    return BBox<dim>(get_vector<dim>(coords, i));
  }
};

namespace sample {  // namespace helps with name resolution in reduction identity 
   template< int N >
   struct bboxWrap {
     BBox<N> box;
  
     KOKKOS_INLINE_FUNCTION   // Default constructor - Initialize to 0's
     bboxWrap() {
       const auto zero = zero_vector<N>();
       box.min = zero;
       box.max = zero;
     }
     KOKKOS_INLINE_FUNCTION   // Copy Constructor
     bboxWrap(const bboxWrap & rhs) { 
       box = rhs.box;
     }
     KOKKOS_INLINE_FUNCTION   // add operator
     bboxWrap& operator += (const bboxWrap& src) {
       box = unite(box,src.box);
       return *this;
     } 
     KOKKOS_INLINE_FUNCTION   // volatile add operator 
     void operator += (const volatile bboxWrap& src) volatile {
       box = unite(box,src.box);
     }
   };
   typedef bboxWrap<3> BB3;  // used to simplify code below
}
namespace Kokkos { //reduction identity must be defined in Kokkos namespace
   template<>
   struct reduction_identity< sample::BB3 > {
      KOKKOS_FORCEINLINE_FUNCTION static sample::BB3 sum() {
         return sample::BB3();
      }
   };
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  const int npts = 3;
  Reals init = {0,0,0,1,2,3,-1,-2,-3};
  Vector<3> min = {-1,-2,-3};
  Vector<3> max = {1,2,3};
  BBox<3> gold(min,max);

  auto transform = GetBBoxOp2<3>(init);

  sample::BB3 tr;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<>(0, npts),
    KOKKOS_LAMBDA(const int& i, sample::BB3& update) {
      update.box = transform(i);
  }, Kokkos::Sum<sample::BB3>(tr) );

  std::cout << "min ";
  for(int i=0; i<3; i++)
    std::cout << tr.box.min[i] << " ";
  std::cout << "\n";
  std::cout << "max ";
  for(int i=0; i<3; i++)
    std::cout << tr.box.max[i] << " ";
  std::cout << "\n";

  assert(are_close(tr.box,gold));

  return 0;
}
