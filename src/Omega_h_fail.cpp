#include <csignal>
#include <cstdarg>
#include <iostream>
#include <sstream>

#include <Omega_h_fail.hpp>

#ifdef OMEGA_H_USE_DWARF
#define BACKWARD_HAS_DWARF 1
#endif
#include <backward/backward.hpp>
namespace backward {
cfile_streambuf::int_type cfile_streambuf::underflow() {
  return traits_type::eof();
}
}  // namespace backward

extern "C" void Omega_h_signal_handler(int s);

namespace Omega_h {

static void print_stacktrace(std::ostream& out, int max_frames) {
  ::backward::StackTrace st;
  st.load_here(std::size_t(max_frames));
  ::backward::Printer p;
  p.print(st, out);
}

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif

void fail(char const* format, ...) {
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
constexpr std::size_t NSIGS = (sizeof(known_signals) / sizeof(known_signals[0]));

void protect() {
  for (std::size_t i = 0; i < NSIGS; ++i) {
    signal(known_signals[i].code, Omega_h_signal_handler);
  }
}

}

extern "C" void Omega_h_signal_handler(int s) {
  static volatile sig_atomic_t already_dying = 0;
  if (already_dying) return;
  already_dying = 1;
  std::stringstream ss;
  for (std::size_t i = 0; i < Omega_h::NSIGS; ++i) {
    if (s == Omega_h::known_signals[i].code) {
      ss << "Omega_h caught signal: " << Omega_h::known_signals[i].name << "\n";
    }
  }
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
