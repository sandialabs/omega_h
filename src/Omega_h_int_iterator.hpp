#ifndef OMEGA_H_INT_ITERATOR_HPP
#define OMEGA_H_INT_ITERATOR_HPP

#include <Omega_h_defines.hpp>
#include <iterator>

namespace Omega_h {

class IntIterator {
  LO i;

 public:
  using value_type = LO;
  using difference_type = LO;
  using reference = LO const&;
  using pointer = LO const*;
  using iterator_category = std::random_access_iterator_tag;
  OMEGA_H_INLINE IntIterator() = default;
  OMEGA_H_INLINE IntIterator(LO i_in) : i(i_in) {}
  OMEGA_H_INLINE bool operator==(IntIterator const& other) const {
    return i == other.i;
  };
  OMEGA_H_INLINE bool operator!=(IntIterator const& other) const {
    return i != other.i;
  };
  OMEGA_H_INLINE reference operator*() const { return i; }
  OMEGA_H_INLINE LO const* operator->() const { return &i; }
  OMEGA_H_INLINE IntIterator& operator++() {
    ++i;
    return *this;
  }
  OMEGA_H_INLINE IntIterator operator++(int) {
    auto ret = *this;
    ++i;
    return ret;
  }
  OMEGA_H_INLINE IntIterator& operator--() {
    --i;
    return *this;
  }
  OMEGA_H_INLINE IntIterator operator--(int) {
    auto ret = *this;
    --i;
    return ret;
  }
  OMEGA_H_INLINE IntIterator& operator+=(difference_type n) {
    i += n;
    return *this;
  }
  OMEGA_H_INLINE IntIterator& operator-=(difference_type n) {
    i -= n;
    return *this;
  }
  OMEGA_H_INLINE IntIterator operator+(difference_type n) const {
    return IntIterator(i + n);
  }
  OMEGA_H_INLINE IntIterator operator-(difference_type n) const {
    return IntIterator(i - n);
  }
  OMEGA_H_INLINE difference_type operator-(IntIterator const& other) const {
    return i - other.i;
  }
  OMEGA_H_INLINE LO operator[](difference_type n) const { return i + n; }
  OMEGA_H_INLINE bool operator<(IntIterator const& other) const {
    return i < other.i;
  }
  OMEGA_H_INLINE bool operator>(IntIterator const& other) const {
    return i > other.i;
  }
  OMEGA_H_INLINE bool operator<=(IntIterator const& other) const {
    return i <= other.i;
  }
  OMEGA_H_INLINE bool operator>=(IntIterator const& other) const {
    return i >= other.i;
  }
};

OMEGA_H_INLINE IntIterator operator+(
    IntIterator::difference_type n, IntIterator it) {
  return it + n;
}

}  // namespace Omega_h

#endif
