#include "Omega_h_int_scan.hpp"

#include "Omega_h_functors.hpp"
#include "Omega_h_scan.hpp"

namespace Omega_h {

template <typename T>
struct ExclScan : public SumFunctor<T> {
  using typename SumFunctor<T>::value_type;
  Read<T> in_;
  Write<LO> out_;
  ExclScan(Read<T> in, Write<LO> out) : in_(in), out_(out) {}
  OMEGA_H_DEVICE void operator()(
      LO i, value_type& update, bool final_pass) const {
    update += in_[i];
    if (final_pass) out_[i + 1] = static_cast<LO>(update);
  }
};

template <typename T>
LOs offset_scan(Read<T> a, std::string const& name) {
  begin_code("offset_scan");
  Write<LO> out(a.size() + 1, name);
  out.set(0, 0);
  parallel_scan(a.size(), ExclScan<T>(a, out), "offset_scan");
  end_code();
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
  parallel_scan(a.size(), FillRight(a), "fill_right");
}

}  // end namespace Omega_h
