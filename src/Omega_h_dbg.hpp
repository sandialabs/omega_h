#ifndef OMEGA_H_DBG_HPP
#define OMEGA_H_DBG_HPP
#include "Omega_h_comm.hpp"
#include "Omega_h_fail.hpp"
#include <iostream>
#include <sstream>

#define OMEGA_H_DBG 1

extern Omega_h::Comm *DBG_COMM;

#define ERROUTFL(a) do { std::ostringstream oss; oss << __FILE__ << ":" << __LINE__ << " :dbg: " << a; OMEGA_H_CHECK(false); } while(0)

#if OMEGA_H_DBG

#define TASK_0_cout if(DBG_COMM && 0 == DBG_COMM->rank()) std::cout
//#define TASK_0_cout std::cout

#define TRACK0(a) do { std::ostringstream oss; oss << __FILE__ << ":" << __LINE__ << " :dbg: " << #a << ": " << a; TASK_0_cout << oss.str() << std::endl; } while(0)
#define TRACK1(a) do { std::ostringstream oss; oss << __FILE__ << ":" << __LINE__ << " :dbg: " << a; TASK_0_cout << oss.str() << std::endl; } while(0)
// FIXME
#define TRACK() TRACK0("srk")

#define TRACK2(a) do { std::ostringstream oss; oss << __FILE__ << ":" << __LINE__ << " :dbg: " << #a << ": " << a; \
                        TASK_0_cout << oss.str() << std::endl; if(DBG_COMM) DBG_COMM->barrier(); } while(0)

#define PCOUTFL(a) do {                                                 \
    if (DBG_COMM) {                                                     \
      for (int irank=0; irank < DBG_COMM->size(); irank++) {            \
        if (DBG_COMM->rank() == irank) {                                \
          std::cout << "P" << irank << ": " << __FILE__ << ":" << __LINE__ << " :dbg: " << a; \
        }                                                               \
        DBG_COMM->barrier();                                            \
      }                                                                 \
    }                                                                   \
  } while(0)

#define PXOUTFL(a) do {                                                 \
    if (DBG_COMM) {                                                     \
      std::cout << "P" << DBG_COMM->rank() << ": " << __FILE__ << ":" << __LINE__ << " :dbg: " << a; \
    }                                                                   \
  } while(0)

namespace Omega_h {

#define DO_STACKTRACE_POS 1
// for more efficiency, less info, use the following:
#define DO_STACKTRACE_POS_USE_INT 0

#if DO_STACKTRACE_POS
#  if DO_STACKTRACE_POS_USE_INT
#    define STACKTRACE_POS_I(i) do { Stacktrace::s_position = (i) + __LINE__; } while(0)
#    define STACKTRACE_POS() STACKTRACE_POS_I(0)
#  else
#    define STACKTRACE_POS() do { std::ostringstream str; str << __FILE__ << ":" << __LINE__; Stacktrace::s_position = str.str(); } while(0)
#  endif

#else
#  define STACKTRACE_POS() do { } while(0)
#endif

class Stacktrace {
public:
#if DO_STACKTRACE_POS
#if DO_STACKTRACE_POS_USE_INT
  static int s_position;
#else
  static std::string s_position;
#endif
#endif
  static void print_stacktrace(size_t sz=50, const std::string& msg="");
  static std::string stacktrace(size_t sz=50, const std::string& msg="");
  static std::string demangle(const std::string& x);
  static std::string demangled_stacktrace(size_t sz=50, bool also_mangled=false, const std::string& msg="");
};

std::string stacktrace(bool print = false);
}

#else // OMEGA_H_DBG

#define TASK_0_cout std::cout

#define TRACK0(a) do {  } while(0)
#define TRACK1(a) do {  } while(0)
// FIXME
#define TRACK() TRACK0("srk")

#define TRACK2(a) do {  } while(0)

#define PCOUTFL(a) do { } while(0)
#define PXOUTFL(a) do { } while(0)

#endif // OMEGA_H_DBG

#endif
