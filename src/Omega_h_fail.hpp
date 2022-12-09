#ifndef OMEGA_H_FAIL_HPP
#define OMEGA_H_FAIL_HPP

#include <Omega_h_config.h>
#ifdef OMEGA_H_ENABLE_DEMANGLED_STACKTRACE
#include <Omega_h_stacktrace.hpp>
#endif
#include <cassert>
#include <iostream>
#include <cstdio>
#include <csignal>
#include <cstdarg>

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
void fail(char const* format, ...);

}  // namespace Omega_h

#define Omega_h_fail Omega_h::fail

#if defined(OMEGA_H_USE_CUDA) && (defined(__clang__) || defined(_MSC_VER))
#if defined(NDEBUG)
#  define OMEGA_H_CHECK(cond) static_cast<void>(cond)
#else
#  define OMEGA_H_CHECK(cond) assert(static_cast<int>(static_cast<bool>(cond)))
#endif
#elif defined(__CUDA_ARCH__)
#if defined(NDEBUG)
#  define OMEGA_H_CHECK(cond) static_cast<void>(cond)
#else
#  define OMEGA_H_CHECK(cond) assert(static_cast<int>(static_cast<bool>(cond)))
#elif defined(OMEGA_H_USE_KOKKOS) && \
      defined(SYCL_LANGUAGE_VERSION) && \
      defined (__INTEL_LLVM_COMPILER)
#define OMEGA_H_CHECK(cond) assert(cond)
#endif
#else
#    define OMEGA_H_CHECK(cond)                                                    \
      ((cond) ? ((void)0)                                                          \
       : Omega_h::fail("assertion %s failed at %s +%d\n", #cond, __FILE__, __LINE__))
#endif

#ifndef NDEBUG
#  define OMEGA_H_CHECK_MSG(cond, a) do {                   \
    if (!(cond)) std::cout << "ERROR: " << a << std::endl;  \
    OMEGA_H_CHECK(cond) ;                                   \
  } while(0)
#  define OMEGA_H_CHECK_PRINTF(cond, format, ...) do {  \
    if (!(cond)) {                                      \
      std::printf("ERROR: " format "\n", __VA_ARGS__);  \
    }                                                   \
    OMEGA_H_CHECK(cond) ;                               \
  } while(0)
#  define OMEGA_H_CHECK_OP(l,op,r) do {                                 \
    if (!((l) op (r))) {                                                \
      std::printf("ERROR: !( %s %s %s ) : ! ( %s %s %s )\n", #l, #op, #r, std::to_string(l).c_str(), #op, std::to_string(r).c_str()); \
    }                                                                   \
    OMEGA_H_CHECK((l) op (r));                                          \
  } while(0)
#  ifdef OMEGA_H_ENABLE_DEMANGLED_STACKTRACE
#    define OMEGA_H_CHECK_OP_1(l,op,r)                                  \
  ((!((l) op (r))) ? ((void)0)                                          \
   : Omega_h::fail("assertion %s %s %s [%d %s %d] failed at %s +%d\n%s\n", #l, #op, #r, l, #op, r, __FILE__, __LINE__))
#  endif

#else // NDEBUG
#  define OMEGA_H_CHECK_MSG(cond, a) OMEGA_H_CHECK(cond)
#  define OMEGA_H_CHECK_PRINTF(cond, format, ...) OMEGA_H_CHECK(cond)
#  define OMEGA_H_CHECK_OP(l,op,r) OMEGA_H_CHECK((l) op (r))
#endif // NDEBUG

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
