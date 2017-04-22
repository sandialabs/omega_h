#include "Omega_h_tag.hpp"

#include "Omega_h_internal.hpp"

#include <regex>

namespace Omega_h {

/* basically just checks for a valid C identifier */
void check_tag_name(std::string const& name) {
  auto id_regex_def = "[_a-zA-Z][_a-zA-Z0-9]*";
  auto id_regex = std::regex(id_regex_def);
  if (!std::regex_match(name, id_regex)) {
    Omega_h_fail("\"%s\" is not a valid tag name.\n"
        "tags must match this regex: %s\n",
        name.c_str(), id_regex_def);
  }
}

TagBase::TagBase(std::string const& name, Int ncomps)
    : name_(name), ncomps_(ncomps) {
  check_tag_name(name);
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
  CHECK(is<T>(t));
  return dynamic_cast<Tag<T> const*>(t);
}

template <typename T>
Tag<T>* as(TagBase* t) {
  CHECK(is<T>(t));
  return dynamic_cast<Tag<T>*>(t);
}

template <typename T>
Tag<T>::Tag(std::string const& name, Int ncomps) : TagBase(name, ncomps) {}

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
  template Tag<T> const* as<T>(TagBase const* t);                              \
  template Tag<T>* as<T>(TagBase * t);                                         \
  template class Tag<T>;
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
