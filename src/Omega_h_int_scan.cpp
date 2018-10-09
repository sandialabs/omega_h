#include "Omega_h_int_scan.hpp"

#include "Omega_h_functors.hpp"
#include "Omega_h_scan.hpp"

namespace Omega_h {

template <typename T>
LOs offset_scan(Read<T> a, std::string const& name) {
  OMEGA_H_TIME_FUNCTION;
  Write<LO> out(a.size() + 1, name);
  out.set(0, 0);
  auto const first = a.begin();
  auto const last = a.end();
  auto const result = out.begin() + 1;
  auto const init = LO(0);
  auto const op = plus<LO>();
  auto transform = identity<LO>();
  transform_inclusive_scan(first, last, result, init, op, std::move(transform));
  return out;
}

template LOs offset_scan(Read<I8> a, std::string const& name);
template LOs offset_scan(Read<I32> a, std::string const& name);

struct FillRight : public MaxFunctor<LO> {
  Write<LO> a_;
  FillRight(Write<LO> a) : a_(a) {}
  OMEGA_H_DEVICE void operator()(
      LO i, value_type& update, bool final_pass) const {
    if (a_[i] > update) update = a_[i];
    if (final_pass && (a_[i] == -1)) a_[i] = static_cast<LO>(update);
  }
};

void fill_right(Write<LO> a) {
  OMEGA_H_TIME_FUNCTION;
  parallel_scan(a.size(), FillRight(a), "fill_right");
}

}  // end namespace Omega_h
