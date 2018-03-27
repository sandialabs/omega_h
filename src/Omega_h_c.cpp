#include <csignal>
#include <cstdarg>
#include <iostream>
#include <sstream>

#include "Omega_h_c.h"
#include "Omega_h_control.hpp"

extern "C" {

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif

void Omega_h_fail(char const* format, ...) {
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  fprintf(stderr, "\n");
  abort();
}

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic pop
#endif

static struct {
  int code;
  const char* name;
} const known_signals[] = {{SIGSYS, "bad system call"},
    {SIGTSTP, "terminal stop"}, {SIGQUIT, "quit"}, {SIGHUP, "hangup"},
    {SIGABRT, "abort"}, {SIGTERM, "termination"},
    {SIGSEGV, "segmentation fault"}, {SIGINT, "interrupt"},
    {SIGILL, "illegal instruction"}, {SIGFPE, "floating point exception"}};
constexpr auto NSIGS = (sizeof(known_signals) / sizeof(known_signals[0]));

void Omega_h_signal_handler(int s) {
  static volatile sig_atomic_t already_dying = 0;
  if (already_dying) return;
  already_dying = 1;
  std::stringstream ss;
  for (size_t i = 0; i < NSIGS; ++i)
    if (s == known_signals[i].code)
      ss << "Omega_h caught signal: " << known_signals[i].name << "\n";
  Omega_h::print_stacktrace(ss, 64);
  auto str = ss.str();
  std::cerr << str;
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
  signal(s, SIG_DFL);
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic pop
#endif
  ::raise(s);
}

void Omega_h_protect() {
  for (size_t i = 0; i < NSIGS; ++i) {
    signal(known_signals[i].code, Omega_h_signal_handler);
  }
}

}  //  extern "C"
