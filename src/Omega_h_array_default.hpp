//
// Created by Matthew McCall on 6/8/23.
//

#ifndef OMEGA_H_ARRAY_DEFAULT_HPP
#define OMEGA_H_ARRAY_DEFAULT_HPP

namespace Omega_h {

template <typename T>
Write<T>::Write() : shared_alloc_() {}

template <typename T>
LO Write<T>::size() const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
  OMEGA_H_CHECK(exists());
#endif
  return static_cast<LO>(shared_alloc_.size() / sizeof(T));
}

template <typename T>
OMEGA_H_DEVICE T& Write<T>::operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
  OMEGA_H_CHECK_OP(0, <=, i);
  OMEGA_H_CHECK_OP(i, <, size());
#endif
  return data()[i];
}

template <typename T>
OMEGA_H_INLINE T* Write<T>::data() const noexcept {
  return static_cast<T*>(shared_alloc_.data());
}

template <typename T>
OMEGA_H_INLINE long Write<T>::use_count() const { return shared_alloc_.alloc->use_count; }

template <typename T>
OMEGA_H_INLINE bool Write<T>::exists() const noexcept {
  return shared_alloc_.data() != nullptr;
}

template <typename T>
inline T const& HostRead<T>::operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
  OMEGA_H_CHECK_OP(0, <=, i);
  OMEGA_H_CHECK_OP(i, <, size());
#endif
#ifdef OMEGA_H_USE_CUDA
  return mirror_[i];
#else
  return read_[i];
#endif
}

template <typename T>
inline T& HostWrite<T>::operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
  OMEGA_H_CHECK_OP(0, <=, i);
  OMEGA_H_CHECK_OP(i, <, size());
#endif
#ifdef OMEGA_H_USE_CUDA
  return mirror_[i];
#else
  return write_[i];
#endif
}

} // namespace Omega_h

#endif  // OMEGA_H_ARRAY_DEFAULT_HPP
