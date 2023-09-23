#ifndef OMEGA_H_MEMORY_HPP
#define OMEGA_H_MEMORY_HPP

#include <map>

namespace Omega_h {

template <typename T>
class SharedRef {
public:
  SharedRef() = default;

  template <typename... Args>
  explicit SharedRef(Args&&... args)
#if defined(OMEGA_H_COMPILING_FOR_HOST)
      : ptr_(new T(std::forward<Args>(args)...)) {
    auto [itr, inserted] = refCount_.insert(std::make_pair(ptr_, 1));
    assert(inserted);
  }
#else
  {
  }
#endif

  OMEGA_H_INLINE SharedRef(const SharedRef& other) {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    if (*this) {
      decrementRefCount();
    }

    if (!other) {
      ptr_ = nullptr;
      return;
    }

    ptr_ = other.ptr_;
    auto itr = refCount_.find(ptr_);
    assert(itr != refCount_.end());
    itr->second++;
#endif
  }

  OMEGA_H_INLINE SharedRef(SharedRef&& other) noexcept {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    if (*this) {
      decrementRefCount();
    }

    if (!other) {
      ptr_ = nullptr;
      return;
    }

    ptr_ = other.ptr_;
    auto itr = refCount_.find(ptr_);
    assert(itr != refCount_.end());
    itr->second++;
#endif
  }

  SharedRef& operator=(const SharedRef& other) {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    if (*this) {
      decrementRefCount();
    }

    if (!other) {
      ptr_ = nullptr;
      return *this;
    }

    ptr_ = other.ptr_;
    auto itr = refCount_.find(ptr_);
    assert(itr != refCount_.end());
    itr->second++;

#endif
    return *this;
  }

  SharedRef& operator=(SharedRef&& other) noexcept {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    if (*this) {
      decrementRefCount();
    }

    if (!other) {
      ptr_ = nullptr;
      return *this;
    }

    ptr_ = other.ptr_;
    auto itr = refCount_.find(ptr_);
    assert(itr != refCount_.end());
    itr->second++;

#endif
    return *this;
  }

  OMEGA_H_INLINE T* get() const {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    return ptr_;
#else
    return nullptr;
#endif
  }

  OMEGA_H_INLINE T* operator->() const {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    return ptr_;
#else
    return nullptr;
#endif
  }

  OMEGA_H_INLINE T& operator*() const {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    return *ptr_;
#else
    return nullptr;
#endif
  }

  OMEGA_H_INLINE ~SharedRef() {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    if (*this) {
      decrementRefCount();
    }
#endif
  }

  OMEGA_H_INLINE explicit operator bool() const {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    return ptr_ != nullptr && (refCount_.find(ptr_) != refCount_.end());
#else
    return false;
#endif
  }

  OMEGA_H_INLINE int use_count() const {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    return *this ? refCount_.find(ptr_)->second : 0;
#else
    return 0;
#endif
  }

private:
  void decrementRefCount() {
#if defined(OMEGA_H_COMPILING_FOR_HOST)
    if (!*this) {
      return;
    }

    auto itr = refCount_.find(ptr_);
    assert(itr != refCount_.end());
    itr->second--;
    if (itr->second == 0) {
      refCount_.erase(itr);
      delete ptr_;
    }
#endif
  }

  T* ptr_ = nullptr;

  static std::map<T*, int> refCount_;
};

template <typename T>
std::map<T*, int> SharedRef<T>::refCount_;

}  // namespace Omega_h

#endif  // OMEGA_H_MEMORY_HPP
