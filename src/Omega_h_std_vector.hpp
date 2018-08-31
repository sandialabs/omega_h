#ifndef OMEGA_H_STD_VECTOR_HPP
#define OMEGA_H_STD_VECTOR_HPP

#include <Omega_h_fail.hpp>
#include <vector>

namespace Omega_h {

/* just some wrappers over std::vector to let us
   do all indexing with int */

template <typename T>
inline int size(std::vector<T> const& v) {
  return int(v.size());
}

template <typename T>
inline typename std::vector<T>::reference at(std::vector<T>& v, int i) {
  OMEGA_H_CHECK(0 <= i);
#if !(defined(__GNUC__) && __GNUC__ < 5)
  OMEGA_H_CHECK(i < int(v.size()));
#endif
  return v[std::size_t(i)];
}

template <typename T>
inline typename std::vector<T>::const_reference at(
    std::vector<T> const& v, int i) {
  OMEGA_H_CHECK(0 <= i);
  OMEGA_H_CHECK(i < int(v.size()));
  return v[std::size_t(i)];
}

template <typename T>
inline void resize(std::vector<T>& v, int n) {
  OMEGA_H_CHECK(0 <= n);
  v.resize(std::size_t(n));
}

template <typename T>
inline void reserve(std::vector<T>& v, int n) {
  OMEGA_H_CHECK(0 <= n);
  v.reserve(std::size_t(n));
}

template <typename T>
inline std::vector<T> make_vector(int n, T const& init_val = T()) {
  return std::vector<T>(std::size_t(n), init_val);
}

}  // namespace Omega_h

#endif
