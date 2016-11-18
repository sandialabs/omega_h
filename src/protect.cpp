#include "protect.hpp"

#include <csignal>
#include <iostream>
#include <sstream>

#include "stacktrace.hpp"

namespace Omega_h {

void osh_signal_handler(int s);

static struct {
  int code;
  const char* name;
} const known_signals[] = {
  {SIGSYS, "bad system call"},
  {SIGTSTP, "terminal stop"},
  {SIGQUIT, "quit"},
  {SIGHUP, "hangup"},
  {SIGABRT, "abort"},
  {SIGTERM, "termination"},
  {SIGSEGV, "segmentation fault"},
  {SIGINT, "interrupt"},
  {SIGILL, "illegal instruction"},
  {SIGFPE, "floating point exception"}
};

#define NSIGS (sizeof(known_signals)/sizeof(known_signals[0]))

void protect() {
  for (size_t i = 0; i < NSIGS; ++i) {
    signal(known_signals[i].code, osh_signal_handler);
  }
}

void osh_signal_handler(int s) {
  static volatile sig_atomic_t already_dying = 0;
  if (already_dying) return;
  already_dying = 1;
  std::stringstream ss;
  for (size_t i = 0; i < NSIGS; ++i)
    if (s == known_signals[i].code)
      ss << "Omega_h caught signal: " << known_signals[i].name << "\n";
  print_stacktrace(ss, 64);
  auto str = ss.str();
  std::cerr << str;
  signal(s, SIG_DFL);
  ::raise(s);
}

#undef NSIGS

}  // end namespace Omega_h
