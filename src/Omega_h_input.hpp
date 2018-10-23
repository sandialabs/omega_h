#ifndef OMEGA_H_INPUT_HPP
#define OMEGA_H_INPUT_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <iosfwd>
#include <Omega_h_defines.hpp>
#include <Omega_h_mesh.hpp>

namespace Omega_h {

struct Input {
  Input* parent;
  bool used;
  virtual ~Input() = default;
  Input(Input const&) = default;
  Input(Input&&) = default;
  Input();
  virtual void out_of_line_virtual_method();
};

std::string get_full_name(Input const& input);

template <class InputType>
bool is_type(Input& input);

template <class InputType>
InputType& as_type(Input& input);

struct InputScalar : public Input {
  std::string str;
  InputScalar(std::string const& str_in);
  InputScalar(InputScalar const&) = default;
  InputScalar(InputScalar&&) = default;
  bool as(std::string& out) const;
  bool as(bool& out) const;
  bool as(double& out) const;
  bool as(int& out) const;
  bool as(long long& out) const;
  template <class T>
  T get() const;
  virtual void out_of_line_virtual_method();
};

struct InputList;

struct InputMap : public Input {
  std::map<std::string, std::shared_ptr<Input>> map;
  InputMap() = default;
  [[noreturn]] InputMap(InputMap const&);
  InputMap(InputMap&&);
  void add(std::string const& name, std::shared_ptr<Input>&& input);
  template <class InputType>
  bool is_input(std::string const& name);
  template <class ScalarType>
  bool is(std::string const& name);
  bool is_map(std::string const& name);
  bool is_list(std::string const& name);
  Input& find_named_input(std::string const& name);
  template <class InputType>
  InputType& use_input(std::string const& name);
  template <class ScalarType>
  ScalarType get(std::string const& name);
  InputMap& get_map(std::string const& name);
  InputList& get_list(std::string const& name);
  template <class ScalarType>
  ScalarType get(std::string const& name, char const* default_value);
  std::string const& name(Input const& input);
  virtual void out_of_line_virtual_method();
};

struct InputList : public Input {
  std::vector<std::shared_ptr<Input>> entries;
  InputList() = default;
  [[noreturn]] InputList(InputList const&);
  InputList(InputList&&);
  void add(std::shared_ptr<Input>&& input);
  LO position(Input const& input);
  LO size();
  Input& at(LO i);
  template <class InputType>
  bool is_input(LO i);
  template <class InputType>
  InputType& use_input(LO i);
  template <class ScalarType>
  bool is(LO i);
  bool is_map(LO i);
  bool is_list(LO i);
  template <class ScalarType>
  ScalarType get(LO i);
  InputMap& get_map(LO i);
  InputList& get_list(LO i);
  virtual void out_of_line_virtual_method();
};

InputMap read_input(std::string const& path);

void update_class_sets(ClassSets* p_sets, InputMap& pl);

void echo_input(std::ostream& stream, Input& input);

void check_unused(Input& input);

#define OMEGA_H_EXPL_INST(InputType) \
extern template bool is_type<InputType>(Input&); \
extern template InputType& as_type<InputType>(Input&); \
extern template bool InputMap::is_input<InputType>(std::string const& name); \
extern template InputType& InputMap::use_input<InputType>(std::string const& name); \
extern template bool InputList::is_input<InputType>(LO i); \
extern template InputType& InputList::use_input(LO i);
OMEGA_H_EXPL_INST(InputScalar)
OMEGA_H_EXPL_INST(InputMap)
OMEGA_H_EXPL_INST(InputList)
#undef OMEGA_H_EXPL_INST

#define OMEGA_H_EXPL_INST(ScalarType) \
extern template ScalarType InputScalar::get<ScalarType>() const; \
extern template bool InputMap::is<ScalarType>(std::string const& name); \
extern template ScalarType InputMap::get<ScalarType>(std::string const& name); \
extern template ScalarType InputMap::get<ScalarType>(std::string const& name, char const* default_value); \
extern template bool InputList::is<ScalarType>(LO i); \
extern template ScalarType InputList::get<ScalarType>(LO i);
OMEGA_H_EXPL_INST(std::string)
OMEGA_H_EXPL_INST(bool)
OMEGA_H_EXPL_INST(double)
OMEGA_H_EXPL_INST(int)
OMEGA_H_EXPL_INST(long long)
#undef OMEGA_H_EXPL_INST

}  // namespace Omega_h

#endif
