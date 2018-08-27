#ifndef OMEGA_H_STRING_HPP
#define OMEGA_H_STRING_HPP

#include <Omega_h_fail.hpp>
#include <string>

namespace Omega_h {

/* some wrappers over std::string to let us
   do all indexing with int */

inline int size(std::string const& s) { return int(s.size()); }

inline typename std::string::reference at(std::string& s, int i) {
  OMEGA_H_CHECK(0 <= i);
  OMEGA_H_CHECK(i < int(s.size()));
  return s[std::size_t(i)];
}

inline typename std::string::const_reference at(std::string const& s, int i) {
  OMEGA_H_CHECK(0 <= i);
  OMEGA_H_CHECK(i < int(s.size()));
  return s[std::size_t(i)];
}

}  // namespace Omega_h

#endif
