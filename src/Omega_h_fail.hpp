#ifndef OMEGA_H_FAIL_H
#define OMEGA_H_FAIL_H

#include <cassert>
#include <Omega_h_config.h>

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
void fail(char const* format, ...)
    __attribute__((noreturn, format(printf, 1, 2)));

}

#define Omega_h_fail Omega_h::fail

#if defined(OMEGA_H_USE_CUDA) && defined(__clang__)
#define OMEGA_H_CHECK(cond) assert(cond)
#elif defined(__CUDA_ARCH__)
#define OMEGA_H_CHECK(cond) assert(cond)
#else
#define OMEGA_H_CHECK(cond)                                                    \
  ((cond) ? ((void)0)                                                          \
          : Omega_h::fail(                                                      \
                "assertion %s failed at %s +%d\n", #cond, __FILE__, __LINE__))
#endif

#if defined( __clang__ ) && !defined(OMEGA_H_USE_CUDA)
#define OMEGA_H_NORETURN(x) OMEGA_H_CHECK(false)
#else
#define OMEGA_H_NORETURN(x)                                                    \
  do {                                                                         \
    OMEGA_H_CHECK(false);                                                             \
    return x;                                                                  \
  } while (false)
#endif

#endif
