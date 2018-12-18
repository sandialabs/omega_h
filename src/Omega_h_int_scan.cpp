#include "Omega_h_int_scan.hpp"

#include "Omega_h_cast_iterator.hpp"
#include "Omega_h_functors.hpp"
#include "Omega_h_profile.hpp"
#include "Omega_h_scan.hpp"

namespace Omega_h {

template <typename T>
LOs offset_scan(Read<T> a, std::string const& name) {
  OMEGA_H_TIME_FUNCTION;
  Write<LO> out(a.size() + 1, name);
  out.set(0, 0);
  auto const first = CastIterator<LO, T>(a.begin());
  auto const last = CastIterator<LO, T>(a.end());
  auto const result = out.begin() + 1;
  inclusive_scan(first, last, result);
  return out;
}

template LOs offset_scan(Read<I8> a, std::string const& name);
template LOs offset_scan(Read<I32> a, std::string const& name);

void fill_right(Write<LO> a) {
  OMEGA_H_TIME_FUNCTION;
  auto const first = a.begin();
  auto const last = a.end();
  auto const result = a.begin();
  auto const op = maximum<LO>();
  auto transform = identity<LO>();
  transform_inclusive_scan(first, last, result, op, std::move(transform));
}

}  // end namespace Omega_h
