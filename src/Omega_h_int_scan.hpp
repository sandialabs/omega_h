#ifndef OMEGA_H_INT_SCAN_HPP
#define OMEGA_H_INT_SCAN_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

template <typename T>
LOs offset_scan(Read<T> a, std::string const& name = "");

extern template LOs offset_scan(Read<I8> a, std::string const& name);
extern template LOs offset_scan(Read<I32> a, std::string const& name);

/* given an array whose values are
   either non-negative or (-1), and whose
   non-negative values are sorted in
   increasing order, replaces all (-1)
   entries with the nearest non-negative
   value to the left */
void fill_right(Write<LO> a);

}  // end namespace Omega_h

#endif
