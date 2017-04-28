#ifndef OMEGA_H_SORT_HPP
#define OMEGA_H_SORT_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

template <typename T>
LOs sort_by_keys(Read<T> keys, Int width = 1);

#define OMEGA_H_INST_DECL(T) extern template LOs sort_by_keys(Read<T> keys, Int width);
OMEGA_H_INST_DECL(LO)
OMEGA_H_INST_DECL(GO)
#undef OMEGA_H_INST_DECL

}  // end namespace Omega_h

#endif
