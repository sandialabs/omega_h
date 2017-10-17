#ifndef OMEGA_H_FUNCTORS_HPP
#define OMEGA_H_FUNCTORS_HPP

#include <Omega_h_scalar.hpp>

namespace Omega_h {

/* Kokkos requires reduction value types to be at least sizeof(int).
   This class encapsulates the choice to use a wider integer if needed.
 */
template <typename T>
struct StandinTraits {
  typedef T type;
};

template <>
struct StandinTraits<I8> {
  typedef I32 type;
};

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
  typedef typename StandinTraits<T>::type value_type;
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
  typedef typename StandinTraits<T>::type value_type;
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
  using value_type = typename StandinTraits<T>::type;
  using input_type = T;
  OMEGA_H_INLINE void init(value_type& update) const { update = 0; }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = update + input;
  }
};

}  // end namespace Omega_h

#endif
