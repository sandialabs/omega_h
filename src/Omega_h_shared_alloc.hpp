#ifndef OMEGA_H_SHARED_ALLOC_HPP
#define OMEGA_H_SHARED_ALLOC_HPP

#include <Omega_h_macros.h>
#include <cstddef>
#include <string>
#include <vector>

// DEBUG!
#include <iostream>

namespace Omega_h {

class Library;

struct Allocs;

extern bool entering_parallel;
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
  ~Alloc();
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
  void* alloc_ptr;
  void* direct_ptr;
  OMEGA_H_INLINE SharedAlloc() noexcept : alloc_ptr(nullptr), direct_ptr(nullptr) {
  }
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
  OMEGA_H_INLINE Alloc* get_alloc() const noexcept { return static_cast<Alloc*>(alloc_ptr); }
  OMEGA_H_INLINE void copy(SharedAlloc const& other) noexcept {
    alloc_ptr = other.alloc_ptr;
#ifndef __CUDA_ARCH__
    if (alloc_ptr && (!(reinterpret_cast<std::uintptr_t>(alloc_ptr) & FREE_MASK))) {
      // allocated
      if (entering_parallel) {
        std::cerr << "entering parallel!\n";
        std::cerr << "get_alloc()->size " << get_alloc()->size << '\n';
        std::cerr << "std::uintptr_t(get_alloc()->size) " << std::uintptr_t(get_alloc()->size) << '\n';
        std::cerr << "std::uintptr_t(get_alloc()->size) << 3 " << (std::uintptr_t(get_alloc()->size) << 3) << '\n';
        std::cerr << "(std::uintptr_t(get_alloc()->size) << 3) | IN_PARALLEL " << ((std::uintptr_t(get_alloc()->size) << 3) | IN_PARALLEL) << '\n';
        std::cerr << "reinterpret_cast<void*>(std::uintptr_t(get_alloc()->size) << 3) | IN_PARALLEL " << reinterpret_cast<void*>((std::uintptr_t(get_alloc()->size) << 3) | IN_PARALLEL) << '\n';
        alloc_ptr = reinterpret_cast<void*>(
            (std::uintptr_t(get_alloc()->size) << 3) | IN_PARALLEL);
        std::cerr << "alloc_ptr " << alloc_ptr << '\n';
      } else {
        ++(get_alloc()->use_count);
      }
    }
#endif
    direct_ptr = other.direct_ptr;
  }
  OMEGA_H_INLINE SharedAlloc(SharedAlloc const& other) {
    copy(other);
  }
  OMEGA_H_INLINE void move(SharedAlloc&& other) {
    alloc_ptr = other.alloc_ptr;
    direct_ptr = other.direct_ptr;
#ifndef __CUDA_ARCH__
    if (alloc_ptr && (!(reinterpret_cast<std::uintptr_t>(alloc_ptr) & FREE_MASK))) {
      // allocated
      if (entering_parallel) {
        --(get_alloc()->use_count);
        alloc_ptr = reinterpret_cast<Alloc*>(
            (std::uintptr_t(get_alloc()->size) << 3) | IN_PARALLEL);
      }
    }
#endif
    other.alloc_ptr = nullptr;
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
    if (alloc_ptr && (!(reinterpret_cast<std::uintptr_t>(alloc_ptr) & FREE_MASK))) {
      // allocated
      --(get_alloc()->use_count);
      if (get_alloc()->use_count == 0) delete get_alloc();
    }
#endif
  }
  OMEGA_H_INLINE std::size_t size() const {
#ifndef __CUDA_ARCH__
    if (!(reinterpret_cast<std::uintptr_t>(alloc_ptr) & IN_PARALLEL)) {
#if defined (__GNUC__) && (__GNUC__ >= 7) && (!defined (__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#endif
      return get_alloc()->size;
#if defined (__GNUC__) && (__GNUC__ >= 7) && (!defined (__clang__))
#pragma GCC diagnostic pop
#endif
    }
#endif
    return reinterpret_cast<std::uintptr_t>(alloc_ptr) >> 3;
  }
  OMEGA_H_INLINE int maybe_identity_index(int i) {
    if (reinterpret_cast<std::uintptr_t>(alloc_ptr) == IS_IDENTITY) {
      return i;
    }
    return static_cast<int*>(direct_ptr)[i];
  }
  OMEGA_H_INLINE void* data() const { return direct_ptr; }
  static SharedAlloc identity(std::size_t size_in);
};

}  // namespace Omega_h

#endif
