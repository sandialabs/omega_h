#ifndef OMEGA_H_TAG_HPP
#define OMEGA_H_TAG_HPP

#include <Omega_h_array.hpp>
#include <string>

namespace Omega_h {

class TagBase {
 public:
  TagBase(std::string const& name, Int ncomps, Int xfer, Int outflags);
  virtual ~TagBase();
  std::string const& name() const;
  Int ncomps() const;
  Int xfer() const;
  Int outflags() const;
  virtual Omega_h_Type type() const = 0;

 private:
  std::string name_;
  Int ncomps_;
  Int xfer_;
  Int outflags_;
};

template <typename T>
class Tag : public TagBase {
 public:
  Tag(std::string const& name, Int ncomps, Int xfer, Int outflags);
  Read<T> array() const;
  void set_array(Read<T> array);
  virtual Omega_h_Type type() const override;

 private:
  Read<T> array_;
};

template <typename T>
bool is(TagBase const* t);

template <typename T>
Tag<T> const* to(TagBase const* t);
template <typename T>
Tag<T>* to(TagBase* t);

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template Read<T> Tag<T>::array() const;                               \
  extern template bool is<T>(TagBase const* t);                                \
  extern template Tag<T> const* to<T>(TagBase const* t);                       \
  extern template Tag<T>* to<T>(TagBase * t);                                  \
  extern template class Tag<T>;
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

}  // namespace Omega_h

#endif
