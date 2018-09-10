#ifndef OMEGA_H_STD_STACK_HPP
#define OMEGA_H_STD_STACK_HPP

#include <stack>

namespace Omega_h {

template <typename T>
int size(std::stack<T> const& s) {
  return static_cast<int>(s.size());
}

}  // namespace Omega_h

#endif
