#ifndef OMEGA_H_SHARED_ALLOC_HPP
#define OMEGA_H_SHARED_ALLOC_HPP

#include <Omega_h_macros.h>
#include <cstddef>
#include <string>
#include <vector>

namespace Omega_h {

class Library;

struct Allocs;

OMEGA_H_DLL extern bool entering_parallel;
extern Allocs* global_allocs;

void start_tracking_allocations();
void stop_tracking_allocations(Library* lib);

struct Alloc {
  std::size_t size;
  std::string name;
  void* ptr;
  int use_count;
  Alloc* prev;
  Alloc* next;
  Alloc(std::size_t size_in, std::string const& name_in);
  Alloc(std::size_t size_in, std::string&& name_in);
  OMEGA_H_DLL ~Alloc();
  Alloc(Alloc const&) = delete;
  Alloc(Alloc&&) = delete;
  Alloc& operator=(Alloc const&) = delete;
  Alloc& operator=(Alloc&&) = delete;
  void init();
};

struct HighWaterRecord {
  std::string name;
  std::size_t bytes;
};

struct Allocs {
  Alloc* first;
  Alloc* last;
  std::size_t total_bytes;
  std::size_t high_water_bytes;
  std::vector<HighWaterRecord> high_water_records;
};

struct SharedAlloc {
  Alloc* alloc;
  void* direct_ptr;
  OMEGA_H_INLINE SharedAlloc() noexcept : alloc(nullptr), direct_ptr(nullptr) {}
  SharedAlloc(std::size_t size_in, std::string const& name_in);
  SharedAlloc(std::size_t size_in, std::string&& name_in);
  SharedAlloc(std::size_t size_in);
  enum : std::uintptr_t {
    FREE_BIT1 = 0x1,
    FREE_BIT2 = 0x2,
    FREE_BIT3 = 0x4,
    FREE_MASK = 0x7,
    IN_PARALLEL = FREE_BIT1,
    IS_IDENTITY = FREE_BIT2,
  };
  /* possible states:
     1. uninitialized: alloc == nullptr
     2. uninitialized, parallel: alloc = IN_PARALLEL
     3. allocated: alloc == true_alloc
     4. allocated, parallel: alloc == (size << 3) & IN_PARALLEL
     5. identity: alloc = (size << 3) & IS_IDENTITY
     6. identity, parallel: alloc = (size << 3) & IS_IDENTITY & IN_PARALLEL
   */
  OMEGA_H_INLINE void copy(SharedAlloc const& other) noexcept {
    alloc = other.alloc;
#ifndef __CUDA_ARCH__
    if (alloc && (!(reinterpret_cast<std::uintptr_t>(alloc) & FREE_MASK))) {
      // allocated
      if (entering_parallel) {
        alloc = reinterpret_cast<Alloc*>(
            (std::uintptr_t(alloc->size) << 3) | IN_PARALLEL);
      } else {
        ++(alloc->use_count);
      }
    }
#endif
    direct_ptr = other.direct_ptr;
  }
  OMEGA_H_INLINE SharedAlloc(SharedAlloc const& other) { copy(other); }
  OMEGA_H_INLINE void move(SharedAlloc&& other) noexcept {
    alloc = other.alloc;
    direct_ptr = other.direct_ptr;
#ifndef __CUDA_ARCH__
    if (alloc && (!(reinterpret_cast<std::uintptr_t>(alloc) & FREE_MASK))) {
      // allocated
      if (entering_parallel) {
        --(alloc->use_count);
        alloc = reinterpret_cast<Alloc*>(
            (std::uintptr_t(alloc->size) << 3) | IN_PARALLEL);
      }
    }
#endif
    other.alloc = nullptr;
    other.direct_ptr = nullptr;
  }
  OMEGA_H_INLINE SharedAlloc(SharedAlloc&& other) { move(std::move(other)); }
  OMEGA_H_INLINE SharedAlloc& operator=(SharedAlloc const& other) {
    clear();
    copy(other);
    return *this;
  }
  OMEGA_H_INLINE SharedAlloc& operator=(SharedAlloc&& other) {
    clear();
    move(std::move(other));
    return *this;
  }
  OMEGA_H_INLINE ~SharedAlloc() { clear(); }
  OMEGA_H_INLINE void clear() {
#ifndef __CUDA_ARCH__
    if (alloc && (!(reinterpret_cast<std::uintptr_t>(alloc) & FREE_MASK))) {
      // allocated
      --(alloc->use_count);
      if (alloc->use_count == 0) delete alloc;
    }
#endif
  }
  OMEGA_H_INLINE std::size_t size() const noexcept {
#ifndef __CUDA_ARCH__
    if (!(reinterpret_cast<std::uintptr_t>(alloc) & IN_PARALLEL)) {
#if defined(__GNUC__) && (__GNUC__ >= 7) && (!defined(__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#endif
      return alloc->size;
#if defined(__GNUC__) && (__GNUC__ >= 7) && (!defined(__clang__))
#pragma GCC diagnostic pop
#endif
    }
#endif
    return reinterpret_cast<std::uintptr_t>(alloc) >> 3;
  }
  OMEGA_H_INLINE int maybe_identity_index(int i) const noexcept {
    if (reinterpret_cast<std::uintptr_t>(alloc) == IS_IDENTITY) {
      return i;
    }
    return static_cast<int*>(direct_ptr)[i];
  }
  OMEGA_H_INLINE void* data() const noexcept { return direct_ptr; }
  static SharedAlloc identity(std::size_t size_in);
};

}  // namespace Omega_h

#endif
