#include "tag.hpp"

namespace Omega_h {

TagBase::TagBase(std::string const& name, Int ncomps, Int xfer)
    : name_(name), ncomps_(ncomps), xfer_(xfer) {}

TagBase::~TagBase() = default;

std::string const& TagBase::name() const { return name_; }

Int TagBase::ncomps() const { return ncomps_; }

Int TagBase::xfer() const { return xfer_; }

template <typename T>
bool is(TagBase const* t) {
  return nullptr != dynamic_cast<Tag<T> const*>(t);
}

template <typename T>
Tag<T> const* to(TagBase const* t) {
  CHECK(is<T>(t));
  return dynamic_cast<Tag<T> const*>(t);
}

template <typename T>
Tag<T>* to(TagBase* t) {
  CHECK(is<T>(t));
  return dynamic_cast<Tag<T>*>(t);
}

template <typename T>
Tag<T>::Tag(std::string const& name, Int ncomps, Int xfer)
    : TagBase(name, ncomps, xfer) {}

template <typename T>
Read<T> Tag<T>::array() const {
  return array_;
}

template <typename T>
void Tag<T>::set_array(Read<T> array) {
  array_ = array;
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
  template Tag<T> const* to<T>(TagBase const* t);                              \
  template Tag<T>* to<T>(TagBase * t);                                         \
  template class Tag<T>;
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
