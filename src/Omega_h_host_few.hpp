#ifndef HOST_FEW_HPP
#define HOST_FEW_HPP

#include "Omega_h_internal.hpp"

namespace Omega_h {

/* this class is different from Omega_h::Few<T,n>
 * because it avoids dealing with volatile
 * operations, thus saving its members from having
 * to implement them. */
template <typename T, Int n>
class HostFew {
  T array_[n];

 public:
  enum { size = n };
  INLINE T& operator[](Int i) { return array_[i]; }
  INLINE T const& operator[](Int i) const { return array_[i]; }
  INLINE HostFew() {}
  HostFew(std::initializer_list<T> l) {
    Int i = 0;
    for (auto v : l) array_[i++] = v;
  }
};

}  // end namespace Omega_h

#endif
