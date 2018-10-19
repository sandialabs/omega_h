#ifndef OMEGA_H_INPUT_HPP
#define OMEGA_H_INPUT_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace Omega_h {

struct Input {
  Input* parent;
  bool used;
  Input();
};

std::string get_full_name(Input& input);

template <class InputType>
bool is_type<InputType>(Input& input);

template <class InputType>
InputType& as_type<InputType>(Input& input);

struct InputScalar : public Input {
  std::string str;
  InputScalar(std::string const& str_in);
  bool as(std::string& out) const;
  bool as(double& out) const;
  bool as(int& out) const;
  bool as(long long& out) const;
  template <class T>
  T get() const;
};

struct InputMap : public Input {
  std::map<std::string, std::unique_ptr<Input>> map;
  void add(std::string const& name, std::unique_ptr<Input>&& input);
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
  ScalarType& get(std::string const& name, char const* default_value);
};

struct InputList : public Input {
  std::vector<std::unique_ptr<Input>> entries;
  void add(std::unique_ptr<Input>&& input);
  LO position(Input& input);
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
};

#define OMEGA_H_EXPL_INST(InputType) \
extern template bool is_type<InputType>(Input&); \
extern template InputType as_type<InputType>(Input&); \
extern template bool InputMap::is_input<InputType>(std::string const& name); \
extern template InputType& InputMap::use_input<InputType>(std::string const& name) \
extern template bool InputList::is_input(LO i); \
extern template InputType& InputList::use_input(LO i);
OMEGA_H_EXPL_INST_INPUT(InputScalar)
OMEGA_H_EXPL_INST_INPUT(InputMap)
OMEGA_H_EXPL_INST_INPUT(InputList)
#undef OMEGA_H_EXPL_INST

#define OMEGA_H_EXPL_INST(ScalarType) \
extern template ScalarType InputScalar::get<ScalarType>() const; \
extern template bool InputMap::is<ScalarType>(std::string const& name); \
extern template ScalarType InputMap::get<ScalarType>(std::string const& name); \
extern template ScalarType& InputMap::get(std::string const& name, char const* default_value); \
extern template bool InputList::is<ScalarType>(LO i); \
extern template ScalarType InputList::get<ScalarType>(LO i);
OMEGA_H_EXPL_INST(std::string)
OMEGA_H_EXPL_INST(double)
OMEGA_H_EXPL_INST(int)
OMEGA_H_EXPL_INST(long long)
#undef OMEGA_H_EXPL_INST

}  // namespace Omega_h

#endif
