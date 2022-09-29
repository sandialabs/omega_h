#include "Omega_h_int_scan.hpp"

#include "Omega_h_cast_iterator.hpp"
#include "Omega_h_functors.hpp"
#include "Omega_h_profile.hpp"
#include "Omega_h_scan.hpp"

namespace Omega_h {

template< typename T>
struct i32Plus {
  KOKKOS_INLINE_FUNCTION
    LO operator()(const T & a, const T & b) const {
      LO al(a);
      LO bl(b);
      LO res = al + bl;
      printf("res al bl %d %d %d\n", res, al, bl);
      return res;
    }
};

template <typename T>
LOs offset_scan(Read<T> a, std::string const& name) {
  OMEGA_H_TIME_FUNCTION;
  Write<LO> out(a.size() + 1, name);
  out.set(0, 0);
#if defined(OMEGA_H_USE_KOKKOS)
  auto outSub = Kokkos::subview(out.view(),std::make_pair(1,out.size()));
  assert(outSub.size()==a.size());
  auto kkOp = i32Plus<T>();
  Kokkos::Experimental::inclusive_scan("omegah_kk_offset_scan", ExecSpace(), a.view(), outSub, kkOp);
  Kokkos::fence();
#else
  auto const first = CastIterator<LO, T>(a.begin());
  auto const last = CastIterator<LO, T>(a.end());
  auto const result = out.begin() + 1;
  inclusive_scan(first, last, result);
#endif
  return out;
}

template LOs offset_scan(Read<I8> a, std::string const& name);
template LOs offset_scan(Read<I32> a, std::string const& name);

void fill_right(Write<LO> a) {
  OMEGA_H_TIME_FUNCTION;
#if defined(OMEGA_H_USE_KOKKOS)
  auto const op = maximum<LO>();
  auto transform = identity<LO>();
  Kokkos::Experimental::transform_inclusive_scan(
    "omegah_kk_fill_right", ExecSpace(),
    Kokkos::Experimental::cbegin(a.view()),
    Kokkos::Experimental::cend(a.view()),
    Kokkos::Experimental::begin(a.view()),
    op, std::move(transform)
  );
#else
  auto const first = a.begin();
  auto const last = a.end();
  auto const result = a.begin();
  auto const op = maximum<LO>();
  auto transform = identity<LO>();
  transform_inclusive_scan(first, last, result, op, std::move(transform));
#endif
}

}  // end namespace Omega_h
