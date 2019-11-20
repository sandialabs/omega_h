#ifndef OMEGA_H_CAST_ITERATOR_HPP
#define OMEGA_H_CAST_ITERATOR_HPP

#include <iterator>

namespace Omega_h {

template <typename To, typename From>
class CastIterator {
  From const* ptr;

 public:
  using value_type = To;
  using difference_type = std::ptrdiff_t;
  using reference = value_type;
  using pointer = value_type const*;
  using iterator_category = std::random_access_iterator_tag;
  inline CastIterator() = default;
  OMEGA_H_INLINE CastIterator(From const* ptr_in) : ptr(ptr_in) {}
  OMEGA_H_INLINE bool operator==(CastIterator const& other) const {
    return ptr == other.ptr;
  }
  OMEGA_H_INLINE bool operator!=(CastIterator const& other) const {
    return ptr != other.ptr;
  }
  OMEGA_H_INLINE reference operator*() const { return *ptr; }
  OMEGA_H_INLINE CastIterator& operator++() {
    ++ptr;
    return *this;
  }
  OMEGA_H_INLINE CastIterator operator++(int) {
    auto const ret = *this;
    ++ptr;
    return ret;
  }
  OMEGA_H_INLINE CastIterator& operator--() {
    --ptr;
    return *this;
  }
  OMEGA_H_INLINE CastIterator operator--(int) {
    auto const ret = *this;
    --ptr;
    return ret;
  }
  OMEGA_H_INLINE CastIterator& operator+=(difference_type n) {
    ptr += n;
    return *this;
  }
  OMEGA_H_INLINE CastIterator& operator-=(difference_type n) {
    ptr -= n;
    return *this;
  }
  OMEGA_H_INLINE CastIterator operator+(difference_type n) const {
    return CastIterator(ptr + n);
  }
  OMEGA_H_INLINE CastIterator operator-(difference_type n) const {
    return CastIterator(ptr - n);
  }
  OMEGA_H_INLINE difference_type operator-(CastIterator const& other) const {
    return ptr - other.ptr;
  }
  OMEGA_H_INLINE value_type operator[](difference_type n) const {
    return *(ptr + n);
  }
  OMEGA_H_INLINE bool operator<(CastIterator const& other) const {
    return ptr < other.ptr;
  }
  OMEGA_H_INLINE bool operator>(CastIterator const& other) const {
    return ptr > other.ptr;
  }
  OMEGA_H_INLINE bool operator<=(CastIterator const& other) const {
    return ptr <= other.ptr;
  }
  OMEGA_H_INLINE bool operator>=(CastIterator const& other) const {
    return ptr >= other.ptr;
  }
};

template <typename To, typename From>
OMEGA_H_INLINE CastIterator<To, From> operator+(
    typename CastIterator<To, From>::difference_type n,
    CastIterator<To, From> it) {
  return it + n;
}

}  // namespace Omega_h

#endif
