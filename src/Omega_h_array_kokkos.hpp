//
// Created by Matthew McCall on 6/8/23.
//

#ifndef OMEGA_H_ARRAY_KOKKOS_HPP
#define OMEGA_H_ARRAY_KOKKOS_HPP

namespace Omega_h {

template <typename T>
OMEGA_H_INLINE Write<T>::Write() : view_() {}

template <typename T>
OMEGA_H_INLINE LO Write<T>::size() const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
  OMEGA_H_CHECK(exists());
#endif
  return static_cast<LO>(view_.size());
}

template <typename T>
OMEGA_H_DEVICE T& Write<T>::operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
  OMEGA_H_CHECK_OP(0, <=, i);
  OMEGA_H_CHECK_OP(i, <, size());
#endif
    return view_(i);
}

template <typename T>
OMEGA_H_INLINE T* Write<T>::data() const noexcept {
    return view_.data();
}

template <typename T>
OMEGA_H_INLINE long Write<T>::use_count() const {
#if !defined(__HIP__) && !defined(__CUDA_ARCH__)
    return manager_.isReferenceCounted() ? manager_.use_count() : view_.use_count();
#else
    return  view_.use_count();
#endif
}

template <typename T>
OMEGA_H_INLINE bool Write<T>::exists() const noexcept {
    return view().data() != nullptr
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE) && (!defined(__CUDA_ARCH__))
           /* deprecated Kokkos behavior: zero-span views have data()==nullptr
            */
           || view().use_count() != 0
#endif
        ;
}

template <typename T>
inline T const& HostRead<T>::operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK_OP(0, <=, i);
    OMEGA_H_CHECK_OP(i, <, size());
#endif
    return mirror_(i);
}

template <typename T>
inline T& HostWrite<T>::operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK_OP(0, <=, i);
    OMEGA_H_CHECK_OP(i, <, size());
#endif
    return mirror_(i);
}

}  // namespace Omega_h

#endif  // OMEGA_H_ARRAY_KOKKOS_HPP
