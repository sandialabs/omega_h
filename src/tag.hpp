#ifndef TAG_HPP
#define TAG_HPP

#include "internal.hpp"

namespace osh {

template <typename T>
bool is(TagBase const* t);

template <typename T>
Tag<T> const* to(TagBase const* t);
template <typename T>
Tag<T>* to(TagBase* t);

#define INST_DECL(T)                                                           \
  extern template bool is<T>(TagBase const* t);                                \
  extern template Tag<T> const* to<T>(TagBase const* t);                       \
  extern template Tag<T>* to<T>(TagBase * t);                                  \
  extern template class Tag<T>;
INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL

}  // end namespace osh

#endif
