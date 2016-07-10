#ifndef TAG_HPP
#define TAG_HPP

namespace osh {

template <typename T>
bool is(TagBase const* t);

template <typename T>
Tag<T> const* to(TagBase const* t);
template <typename T>
Tag<T>* to(TagBase* t);

} //end namespace osh

#endif
