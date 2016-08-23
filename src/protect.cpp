#include "protect.hpp"

#include <csignal>
#include <cstdio>

#include "stacktrace.hpp"

namespace Omega_h {

#define NSIGS 6

void osh_signal_handler(int s);

static struct {
  int code;
  const char* name;
} const known_signals[NSIGS] = {{SIGTERM, "termination"}, {SIGABRT, "abort"},
    {SIGSEGV, "segmentation fault"}, {SIGINT, "interrupt"},
    {SIGILL, "illegal instruction"}, {SIGFPE, "floating point exception"}};

void protect() {
  for (int i = 0; i < NSIGS; ++i) {
    signal(known_signals[i].code, osh_signal_handler);
  }
}

void osh_signal_handler(int s) {
  static volatile sig_atomic_t already_dying = 0;
  if (already_dying) return;
  already_dying = 1;
  for (int i = 0; i < NSIGS; ++i)
    if (s == known_signals[i].code)
      fprintf(stderr, "omega_h caught %s\n", known_signals[i].name);
  print_stacktrace();
  signal(s, SIG_DFL);
  ::raise(s);
}

#undef NSIGS

}  // end namespace Omega_h
