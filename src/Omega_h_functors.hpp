#ifndef OMEGA_H_FUNCTORS_HPP
#define OMEGA_H_FUNCTORS_HPP

#include <Omega_h_scalar.hpp>

namespace Omega_h {

struct AndFunctor {
  typedef I64 value_type;
  OMEGA_H_INLINE void init(value_type& update) const { update = 1; }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = update && input;
  }
};

template <typename T>
struct MaxFunctor {
  typedef promoted_t<T> value_type;
  typedef T input_type;
  OMEGA_H_INLINE void init(value_type& update) const {
    update = ArithTraits<T>::min();
  }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = max2(update, input);
  }
};

template <typename T>
struct MinFunctor {
  typedef promoted_t<T> value_type;
  typedef T input_type;
  OMEGA_H_INLINE void init(value_type& update) const {
    update = ArithTraits<T>::max();
  }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = min2(update, input);
  }
};

template <typename T>
struct SumFunctor {
  using value_type = promoted_t<T>;
  using input_type = T;
  OMEGA_H_INLINE void init(value_type& update) const { update = 0; }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = update + input;
  }
};

}  // end namespace Omega_h

#endif
