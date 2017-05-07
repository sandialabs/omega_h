#include "Omega_h_tag.hpp"

/* regex support is a complete lie prior to GCC 4.9.0 */

#if !defined(__clang__) && defined(__GNUC__)
/* Test for GCC > 3.2.0 */
#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ >= 9))
#define OMEGA_H_COMPILER_SUPPORTS_REGEX
#endif
#else
#define OMEGA_H_COMPILER_SUPPORTS_REGEX
#endif

#ifdef OMEGA_H_COMPILER_SUPPORTS_REGEX
#include <regex>
#endif

namespace Omega_h {

/* basically just checks for a valid C identifier */
void check_tag_name(std::string const& name) {
#ifdef OMEGA_H_COMPILER_SUPPORTS_REGEX
  auto id_regex_def = "[_a-zA-Z][_a-zA-Z0-9]*";
  auto id_regex = std::regex(id_regex_def);
  if (!std::regex_match(name, id_regex)) {
    Omega_h_fail(
        "\"%s\" is not a valid tag name.\n"
        "tags must match this regex: %s\n",
        name.c_str(), id_regex_def);
  }
#else
  (void) name;
#endif
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
  OMEGA_H_CHECK(is<T>(t));
  return dynamic_cast<Tag<T> const*>(t);
}

template <typename T>
Tag<T>* as(TagBase* t) {
  OMEGA_H_CHECK(is<T>(t));
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
