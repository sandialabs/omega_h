#ifndef OMEGA_H_FAIL_HPP
#define OMEGA_H_FAIL_HPP

#include <Omega_h_config.h>
#include <csignal>
#include <cstdarg>
#include <iostream>
#include <cassert>

#ifdef OMEGA_H_THROW
#include <exception>
#include <string>
#endif

namespace Omega_h {

#ifdef OMEGA_H_THROW
struct exception : public std::exception {
  exception(std::string const& msg_in);
  std::string msg;
  const char* what() const noexcept override;
};
#endif

void protect();	
#ifdef _MSC_VER
__declspec(noreturn)
#else
__attribute__((noreturn, format(printf, 1, 2)))
#endif

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif

static inline void fail(char const* format, ...) {
  va_list vlist;
  va_start(vlist, format);
#ifdef OMEGA_H_THROW
  char buffer[2048];
  std::vsnprintf(buffer, sizeof(buffer), format, vlist);
  va_end(vlist);
  throw Omega_h::exception(buffer);
#else
  //std::vfprintf(stderr, format, vlist);
  //va_end(vlist);
  std::abort();
#endif
}

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic pop
#endif

}  // namespace Omega_h

#define Omega_h_fail Omega_h::fail

#if defined(OMEGA_H_USE_CUDA) && defined(__clang__)
#define OMEGA_H_CHECK(cond) assert(cond)
#elif defined(__CUDA_ARCH__)
#define OMEGA_H_CHECK(cond) assert(cond)
#elif defined(OMEGA_H_USE_CUDA) && defined(OMEGA_H_USE_KOKKOS)
#define OMEGA_H_CHECK(cond) assert(cond)
#elif defined(OMPTARGET)
#define OMEGA_H_CHECK(cond) assert(cond)
#elif defined(OMEGA_H_USE_KOKKOS) && \
      defined(SYCL_LANGUAGE_VERSION) && \
      defined (__INTEL_LLVM_COMPILER)
#define OMEGA_H_CHECK(cond) assert(cond)
#else
#define OMEGA_H_CHECK(cond)                                                    \
  ((cond) ? ((void)0)                                                          \
          : Omega_h::fail(                                                     \
                "assertion %s failed at %s +%d\n", #cond, __FILE__, __LINE__))
#endif

#if defined(__clang__) && !defined(OMEGA_H_USE_CUDA)
#define OMEGA_H_NORETURN(x) OMEGA_H_CHECK(false)
#elif defined(OMEGA_H_USE_KOKKOS) && \
      defined(SYCL_LANGUAGE_VERSION) && \
      defined (__INTEL_LLVM_COMPILER)
#define OMEGA_H_NORETURN(x) assert(false)
#else
#define OMEGA_H_NORETURN(x)                                                    \
  do {                                                                         \
    OMEGA_H_CHECK(false);                                                      \
    return x;                                                                  \
  } while (false)
#endif

#endif
