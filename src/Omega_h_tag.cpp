#include "Omega_h_tag.hpp"

namespace Omega_h {

TagBase::TagBase(std::string const& name_in, Int ncomps_in)
    : name_(name_in), ncomps_(ncomps_in) {
  check_tag_name(name_in);
}

TagBase::~TagBase() = default;

std::string const& TagBase::name() const { return name_; }

Int TagBase::ncomps() const { return ncomps_; }

template <typename T>
bool is(TagBase const* t) {
  return nullptr != dynamic_cast<Tag<T> const*>(t);
}

template <typename T>
Tag<T> const* as(TagBase const* t) {
  OMEGA_H_CHECK(is<T>(t));
  return dynamic_cast<Tag<T> const*>(t);
}

template <typename T>
Tag<T>* as(TagBase* t) {
  OMEGA_H_CHECK(is<T>(t));
  return dynamic_cast<Tag<T>*>(t);
}

template <typename T>
Tag<T>::Tag(std::string const& name_in, Int ncomps_in)
    : TagBase(name_in, ncomps_in) {}

template <typename T>
Read<T> Tag<T>::array() const {
  return array_;
}

template <typename T>
void Tag<T>::set_array(Read<T> array_in) {
  array_ = array_in;
}

template <typename T>
struct TagTraits;

template <>
struct TagTraits<I8> {
  static Omega_h_Type type() { return OMEGA_H_I8; }
};

template <>
struct TagTraits<I32> {
  static Omega_h_Type type() { return OMEGA_H_I32; }
};

template <>
struct TagTraits<I64> {
  static Omega_h_Type type() { return OMEGA_H_I64; }
};

template <>
struct TagTraits<Real> {
  static Omega_h_Type type() { return OMEGA_H_F64; }
};

template <typename T>
Omega_h_Type Tag<T>::type() const {
  return TagTraits<T>::type();
}

#define INST(T)                                                                \
  template bool is<T>(TagBase const* t);                                       \
  template Tag<T> const* as<T>(TagBase const* t);                              \
  template Tag<T>* as<T>(TagBase * t);                                         \
  template class Tag<T>;
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
