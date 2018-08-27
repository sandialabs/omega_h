#ifndef OMEGA_H_INPUT_HPP
#define OMEGA_H_INPUT_HPP

#include <string>
#include <memory>
#include <vector>

namespace Omega_h {

struct Input {
  Input():parent(nullptr),used(false) {}
  Input* parent;
  bool used;
};

inline bool is_type<InputType>(Input& input) {
  auto ptr_in = &input;
  auto ptr_out = dynamic_cast<InputType*>(ptr_in);
  return ptr_out != nullptr;
}

inline InputType& as_type<InputType>(Input& input) {
  return dynamic_cast<InputType&>(input);
}

template <class ScalarType>
struct ScalarInput : public Input {
  ScalarInput(ScalarType value_in):value(value_in) {}
  ScalarType value;
};

struct NamedInput : public Input {
  NamedInput(std::string const& name_in, Input* value_in):name(name_in),value(value_in) {}
  std::string name;
  std::unique_ptr<Input> value;
};

struct MapInput : public Input;
struct ListInput : public Input;

template <class InputType>
InputType& get_input(Input& input) {
  auto out = as_type<InputType>(input);
  out.used = true;
  return out;
}

struct MapInput : public Input {
  std::vector<std::unique_ptr<NamedInput>> entries;
  std::rb_tree<std::string, NamedInput*, NameOfNamedInput> by_name;
  void add(NamedInput* named_input) {
    auto it = by_name.upper_bound(named_input->name);
    if (it != by_name.end() && (*it)->name == named_input->name) {
      fail("Tried to add a mapped input value of name \"%s\" that already existed\n", named_input->name.c_str());
    }
    std::unique_ptr<NamedInput> uptr(named_input);
    entries.push_back(std::move(uptr));
    by_name.insert(it, named_input);
  }
  template <class InputType>
  bool is_input(std::string const& name) {
    auto it = by_name.find(name);
    if (it == by_name.end()) return false;
    auto& uptr = *it;
    return is_type(*uptr);
  }
  template <class ScalarType>
  bool is(std::string const& name) {
    return is_input<ScalarInput<ScalarType>>(name);
  }
  bool is_map(std::string const& name) {
    return is_input<MapInput>(name);
  }
  bool is_list(std::string const& name) {
    return is_input<ListInput>(name);
  }
  Input& get_named_input(std::string const& name) {
    auto it = by_name.find(name);
    if (it == by_name.end()) {
      auto s = get_full_name(*this) + name;
      fail("Tried to get named input \"%s\" that doesn't exist\n", s.c_str()); 
    }
    auto& named_uptr = *it;
    auto& value_uptr = named_uptr->value;
    return *value_uptr;
  }
  template <class InputType>
  InputType& get_input(std::string const& name) {
    return Omega_h::get_input<InputType>(get_named_input(name));
  }
  template <class ScalarType>
  ScalarType& get(std::string const& name) {
    return this->get_input<ScalarInput<ScalarType>>(name).value;
  }
  MapInput& get_map(std::string const& name) {
    return this->get_input<MapInput>(name).value;
  }
  ListInput& get_list(std::string const& name) {
    return this->get_input<ListInput>(name).value;
  }
  template <class ScalarType>
  ScalarType& get(std::string const& name, ScalarType const& default_value) {
    if (has_input<ScalarInput<ScalarType>>(name)) return this->get<ScalarType>(name);
    this->add(new NamedInput(name, new ScalarInput<ScalarType>(default_value)));
  }
};

struct ListInput : public Input {
  std::vector<std::unique_ptr<Input>> entries;
  std::size_t position(Input& input) {
    auto it = std::find(entries.begin(), entries.end(), [&](std::unique_ptr<Input> const& uptr)
        { return uptr.get() == &input; });
    OMEGA_H_CHECK(it != entries.end());
    return it - entries.begin();
  }
  LO size() { return LO(entries.size()); }
  template <class InputType>
  bool is_input(LO i) {
    return Omega_h::is_type<InputType>(*(entries[std::size_t(i)]));
  }
  template <class InputType>
  InputType& get_input(LO i) {
    return Omega_h::get_input<InputType>(*(entries[std::size_t(i)]));
  }
  template <class ScalarType>
  ScalarType& get_input(LO i) {
    return Omega_h::get_input<InputType>(*(entries[std::size_t(i)]));
  }
};

std::string get_full_name(Input& input) {
  std::string full_name;
  if (input.parent != nullptr) {
    auto& parent = *(input.parent);
    full_name = get_full_name(parent);
    if (is_type<ListInput>(parent)) {
      auto i = as_type<ListInput>(parent).position(input);
      full_name += "[";
      full_name += std::to_string(i);
      full_name += "]";
    }
  }
  if (is_type<NamedInput>(input)) {
    if (!full_name.empty()) full_name += ".";
    full_name += as_type<NamedInput>(input).name;
  }
}

};

#endif
