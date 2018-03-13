#ifndef OMEGA_H_PRINT_HPP
#define OMEGA_H_PRINT_HPP

#include <Omega_h_array.hpp>
#include <iostream>

namespace Omega_h {

template <class T>
std::ostream& operator<<(std::ostream& stream, HostRead<T> hr); 
template <class T>
std::ostream& operator<<(std::ostream& stream, Read<T> r); 

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template std::ostream& operator<<(std::ostream&, HostRead<T>); \
  extern template std::ostream& operator<<(std::ostream&, Read<T>);
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

}

#endif
