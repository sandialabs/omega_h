#ifndef OMEGA_H_DBG_HPP
#define OMEGA_H_DBG_HPP
#include "Omega_h_comm.hpp"
#include "Omega_h_fail.hpp"
#include <iostream>
#include <sstream>

extern Omega_h::Comm *DBG_COMM;

#define ERROUTFL(a) do { std::ostringstream oss; oss << __FILE__ << ":" << __LINE__ << " :dbg: " << a; OMEGA_H_CHECK(false); } while(0)

#ifdef OMEGA_H_DBG

namespace Omega_h {
inline std::string proc() { std::ostringstream oss; oss << "P" << (DBG_COMM ? DBG_COMM->rank() : 0) << ": "; return oss.str(); }
}

#  define TASK_0_cout if(DBG_COMM && (0 == DBG_COMM->rank())) std::cout

#  define TRACK0(a) do { std::ostringstream oss; oss << __FILE__ << ":" << __LINE__ << " :dbg: " << #a << ": " << a; TASK_0_cout << oss.str() << std::endl; } while(0)
#  define TRACK1(a) do { std::ostringstream oss; oss << __FILE__ << ":" << __LINE__ << " :dbg: " << a; TASK_0_cout << oss.str() << std::endl; } while(0)
#  define TRACK2(a) do { std::ostringstream oss; oss << __FILE__ << ":" << __LINE__ << " :dbg: " << #a << ": " << a; \
                        TASK_0_cout << oss.str() << std::endl; if(DBG_COMM) DBG_COMM->barrier(); } while(0)
#  define TRACK() TRACK0("track")

#  define PCOUTFL(a) do {                                               \
    if (DBG_COMM) {                                                     \
      for (int irank=0; irank < DBG_COMM->size(); irank++) {            \
        if (DBG_COMM->rank() == irank) {                                \
          std::cout << "P" << irank << ": " << __FILE__ << ":" << __LINE__ << " :dbg: " << a; \
        }                                                               \
        DBG_COMM->barrier();                                            \
      }                                                                 \
    }                                                                   \
  } while(0)

#  define PXOUTFL(a) do {                                               \
    if (DBG_COMM) {                                                     \
      std::cout << "P" << DBG_COMM->rank() << ": " << __FILE__ << ":" << __LINE__ << " :dbg: " << a; \
    }                                                                   \
  } while(0)


#else // OMEGA_H_DBG

namespace Omega_h {
inline std::string proc() { return "P0: "; }
}

#  define TASK_0_cout std::cout

#  define TRACK0(a) do { } while(0)
#  define TRACK1(a) do { } while(0)
#  define TRACK() do { } while(0)

#  define TRACK2(a) do { } while(0)

#  define PCOUTFL(a) do { } while(0)
#  define PXOUTFL(a) do { } while(0)

#endif // OMEGA_H_DBG

#endif // OMEGA_H_DBG_HPP
