#include <Omega_h_fail.hpp>
#include <Omega_h_input.hpp>
#include <Omega_h_profile.hpp>
#include <Omega_h_reader.hpp>
#include <Omega_h_yaml.hpp>
#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>

namespace Omega_h {

Input::Input() : parent(nullptr), used(false) {}

void Input::out_of_line_virtual_method() {}

std::string get_full_name(Input const& input) {
  std::string full_name;
  if (input.parent != nullptr) {
    auto& parent = *(input.parent);
    full_name = get_full_name(parent);
    if (is_type<InputList>(parent)) {
      auto i = as_type<InputList>(parent).position(input);
      full_name += "[";
      full_name += std::to_string(i);
      full_name += "]";
    } else if (is_type<InputMap>(parent)) {
      auto& name = as_type<InputMap>(parent).name(input);
      full_name += ".";
      full_name += name;
    }
  }
  return full_name;
}

template <class InputType>
bool is_type(Input& input) {
  auto ptr_in = &input;
  auto ptr_out = dynamic_cast<InputType*>(ptr_in);
  return ptr_out != nullptr;
}

template <class InputType>
InputType& as_type(Input& input) {
  return dynamic_cast<InputType&>(input);
}

InputScalar::InputScalar(std::string const& str_in) : str(str_in) {}

bool InputScalar::as(std::string& out) const {
  out = str;
  return true;
}

bool InputScalar::as(bool& out) const {
  if (str == "true") {
    out = true;
    return true;
  }
  if (str == "false") {
    out = false;
    return true;
  }
  return false;
}

bool InputScalar::as(double& out) const {
  std::istringstream ss(str);
  ss >> std::noskipws >> out;
  return ss.eof() && !ss.fail();
}

bool InputScalar::as(int& out) const {
  std::istringstream ss(str);
  using LL = long long;
  LL val;
  ss >> std::noskipws >> val;
  if (ss.eof() && !ss.fail() && (val >= LL(std::numeric_limits<int>::min())) &&
      (val <= LL(std::numeric_limits<int>::max()))) {
    out = int(val);
    return true;
  }
  return false;
}

bool InputScalar::as(long long& out) const {
  std::istringstream ss(str);
  ss >> std::noskipws >> out;
  return ss.eof() && !ss.fail();
}

template <class T>
T InputScalar::get() const {
  T out;
  if (!as(out)) {
    auto full_name = get_full_name(*this);
    Omega_h_fail(
        "InputScalar \"%s\" string \"%s\" is not interpretable as a %s",
        full_name.c_str(), str.c_str(),
        (std::is_same<T, int>::value
                ? "int"
                : (std::is_same<T, double>::value
                          ? "double"
                          : (std::is_same<T, long long>::value
                                    ? "long long"
                                    : "unknown type"))));
  }
  return out;
}

void InputScalar::out_of_line_virtual_method() {}

InputMapIterator::InputMapIterator(decltype(impl) impl_in) : impl(impl_in) {}

bool InputMapIterator::operator==(InputMapIterator const& other) const
    noexcept {
  return impl == other.impl;
}

bool InputMapIterator::operator!=(InputMapIterator const& other) const
    noexcept {
  return impl != other.impl;
}

InputMapIterator::reference InputMapIterator::operator*() const noexcept {
  return impl->first;
}

InputMapIterator::pointer InputMapIterator::operator->() const noexcept {
  return &(impl->first);
}

InputMapIterator& InputMapIterator::operator++() noexcept {
  ++impl;
  return *this;
}

InputMapIterator InputMapIterator::operator++(int) noexcept { return impl++; }

InputMapIterator& InputMapIterator::operator--() noexcept {
  --impl;
  return *this;
}

InputMapIterator InputMapIterator::operator--(int) noexcept { return impl--; }

InputMap::InputMap(InputMap&& other) : Input(other), map(std::move(other.map)) {
  for (auto& pair : map) (pair.second)->parent = this;
}

InputMap::InputMap(InputMap const& other) : Input(other) {
  Omega_h_fail("InputMap should never actually be copied!\n");
}

void InputMap::add(std::string const& name, std::shared_ptr<Input>&& input) {
  input->parent = this;
  auto const did = map.emplace(name, std::move(input)).second;
  if (!did) {
    fail("tried to add already existing InputMap name \"%s\"\n", name.c_str());
  }
}

template <class InputType>
bool InputMap::is_input(std::string const& name) {
  auto const it = map.find(name);
  if (it == map.end()) return false;
  auto const& sptr = it->second;
  return is_type<InputType>(*sptr);
}

template <class ScalarType>
bool InputMap::is(std::string const& name) {
  auto const it = map.find(name);
  if (it == map.end()) return false;
  auto const& sptr = it->second;
  ScalarType ignored;
  return is_type<InputScalar>(*sptr) && as_type<InputScalar>(*sptr).as(ignored);
}

bool InputMap::is_map(std::string const& name) {
  return this->is_input<InputMap>(name);
}

bool InputMap::is_list(std::string const& name) {
  return this->is_input<InputList>(name);
}

Input& InputMap::find_named_input(std::string const& name) {
  auto it = map.find(name);
  if (it == map.end()) {
    auto s = get_full_name(*this) + "." + name;
    fail("tried to find InputMap entry \"%s\" that doesn't exist\n", s.c_str());
  }
  return *(it->second);
}

template <class InputType>
static InputType& use_input(Input& input) {
  input.used = true;
  return as_type<InputType>(input);
}

template <class InputType>
InputType& InputMap::use_input(std::string const& name) {
  return Omega_h::use_input<InputType>(find_named_input(name));
}

template <class ScalarType>
ScalarType InputMap::get(std::string const& name) {
  return this->use_input<InputScalar>(name).get<ScalarType>();
}

InputMap& InputMap::get_map(std::string const& name) {
  if (!is_map(name)) {
    std::shared_ptr<Input> sptr(new InputMap());
    this->add(name, std::move(sptr));
  }
  return this->use_input<InputMap>(name);
}

InputList& InputMap::get_list(std::string const& name) {
  if (!is_list(name)) {
    std::shared_ptr<Input> sptr(new InputList());
    this->add(name, std::move(sptr));
  }
  return this->use_input<InputList>(name);
}

template <class ScalarType>
ScalarType InputMap::get(std::string const& name, char const* default_value) {
  if (!this->is<ScalarType>(name)) set(name, default_value);
  return this->get<ScalarType>(name);
}

void InputMap::set(std::string const& name, char const* value) {
  std::shared_ptr<Input> sptr(new InputScalar(value));
  this->add(name, std::move(sptr));
}

std::string const& InputMap::name(Input const& input) {
  for (auto& pair : map) {
    if (pair.second.get() == &input) return pair.first;
  }
  OMEGA_H_NORETURN(map.end()->first);
}

void InputMap::out_of_line_virtual_method() {}

InputList::InputList(InputList&& other) : entries(std::move(other.entries)) {
  for (auto& sptr : entries) sptr->parent = this;
}

InputList::InputList(InputList const& other) : Input(other) {
  Omega_h_fail("InputList should never actually be copied!\n");
}

void InputList::add(std::shared_ptr<Input>&& input) {
  input->parent = this;
  entries.push_back(std::move(input));
}

LO InputList::position(Input const& input) {
  auto it = std::find_if(entries.begin(), entries.end(),
      [&](std::shared_ptr<Input> const& sptr) { return sptr.get() == &input; });
  OMEGA_H_CHECK(it != entries.end());
  return LO(it - entries.begin());
}

LO InputList::size() { return LO(entries.size()); }

Input& InputList::at(LO i) { return *(entries.at(std::size_t(i))); }

template <class InputType>
bool InputList::is_input(LO i) {
  return Omega_h::is_type<InputType>(this->at(i));
}

template <class InputType>
InputType& InputList::use_input(LO i) {
  return Omega_h::use_input<InputType>(this->at(i));
}

template <class ScalarType>
bool InputList::is(LO i) {
  Input& input = at(i);
  ScalarType ignored;
  return is_type<InputScalar>(input) && as_type<InputScalar>(input).as(ignored);
}

bool InputList::is_map(LO i) { return this->is_input<InputMap>(i); }

bool InputList::is_list(LO i) { return this->is_input<InputList>(i); }

template <class ScalarType>
ScalarType InputList::get(LO i) {
  return this->use_input<InputScalar>(i).get<ScalarType>();
}

InputMap& InputList::get_map(LO i) { return this->use_input<InputMap>(i); }

InputList& InputList::get_list(LO i) { return this->use_input<InputList>(i); }

void InputList::out_of_line_virtual_method() {}

InputMapIterator InputMap::begin() const noexcept { return map.begin(); }

InputMapIterator InputMap::end() const noexcept { return map.end(); }

void InputMap::remove(std::string const& name) {
  auto const it = map.find(name);
  if (it == map.end())
    Omega_h_fail("InputMap::remove: \"%s\" didn't exist\n", name.c_str());
  map.erase(it);
}

struct NameValue {
  std::string name;
  std::shared_ptr<Input> value;
};

static std::string remove_trailing_whitespace(std::string const& in) {
  std::size_t new_end = 0;
  for (std::size_t ri = 0; ri < in.size(); ++ri) {
    std::size_t i = in.size() - 1 - ri;
    if (in[i] != ' ' && in[i] != '\t') {
      new_end = i + 1;
      break;
    }
  }
  return in.substr(0, new_end);
}

static std::string remove_trailing_whitespace_and_newlines(
    std::string const& in) {
  std::size_t new_end = 0;
  for (std::size_t ri = 0; ri < in.size(); ++ri) {
    std::size_t i = in.size() - 1 - ri;
    if (in[i] != ' ' && in[i] != '\t' && in[i] != '\n' && in[i] != '\r') {
      new_end = i + 1;
      break;
    }
  }
  return in.substr(0, new_end);
}

// http://en.cppreference.com/w/cpp/string/byte/isdigit
static bool my_isdigit(char ch) {
  return std::isdigit(static_cast<unsigned char>(ch));
}

class InputYamlReader : public Reader {
 public:
  InputYamlReader() : Reader(yaml::ask_reader_tables()) {}
  ~InputYamlReader() override;

 protected:
  enum { TRIM_NORMAL, TRIM_DASH };
  any at_shift(int token, std::string& text) override final {
    switch (token) {
      case yaml::TOK_NEWLINE: {
        return any(std::move(text));
      }
      case yaml::TOK_SPACE:
      case yaml::TOK_OTHER: {
        return text.at(0);
      }
    }
    return any();
  }
  any at_reduce(int prod, std::vector<any>& rhs) override final {
    switch (prod) {
      case yaml::PROD_DOC:
      case yaml::PROD_DOC2: {
        std::size_t offset = prod == yaml::PROD_DOC2 ? 1 : 0;
        OMEGA_H_CHECK(!rhs.at(offset).empty());
        auto& result_any = rhs.at(offset);
        OMEGA_H_CHECK(result_any.type() == typeid(InputMap));
        auto& map = any_cast<InputMap&>(result_any);
        map.used = true;
        return any(std::move(result_any));
      }
      case yaml::PROD_TOP_BMAP: {
        return any(std::move(rhs.at(0)));
      }
      case yaml::PROD_TOP_FIRST: {
        if (rhs.at(0).type() == typeid(NameValue)) {
          return map_first_item(rhs.at(0));
        }
        return any();
      }
      case yaml::PROD_TOP_NEXT: {
        if (rhs.at(1).type() == typeid(NameValue)) {
          if (rhs.at(0).type() == typeid(InputMap)) {
            return map_next_item(rhs.at(0), rhs.at(1));
          } else {
            return map_first_item(rhs.at(0));
          }
        } else {
          return any(std::move(rhs.at(0)));
        }
      }
      case yaml::PROD_BMAP_FIRST:
      case yaml::PROD_FMAP_FIRST: {
        OMEGA_H_CHECK(rhs.at(0).type() == typeid(NameValue));
        auto result_any = map_first_item(rhs.at(0));
        OMEGA_H_CHECK(result_any.type() == typeid(InputMap));
        return result_any;
      }
      case yaml::PROD_BMAP_NEXT: {
        return map_next_item(rhs.at(0), rhs.at(1));
      }
      case yaml::PROD_FMAP_NEXT: {
        return map_next_item(rhs.at(0), rhs.at(3));
      }
      case yaml::PROD_BMAP_SCALAR:
      case yaml::PROD_FMAP_SCALAR: {
        return map_item(rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_FMAP_FMAP: {
        return map_item(rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_FMAP_FSEQ: {
        return map_item(rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_BMAP_BSCALAR: {
        return map_item(rhs.at(0), rhs.at(3));
      }
      case yaml::PROD_BMAP_BVALUE: {
        return map_item(rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_BVALUE_EMPTY: {
        return InputMap();
      }
      case yaml::PROD_BVALUE_BMAP: {
        return any(std::move(rhs.at(1)));
      }
      case yaml::PROD_BVALUE_BSEQ: {
        return any(std::move(rhs.at(1)));
      }
      case yaml::PROD_BMAP_FMAP: {
        return map_item(rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_BMAP_FSEQ: {
        return map_item(rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_BSEQ_FIRST: {
        return seq_first_item(rhs.at(0));
      }
      case yaml::PROD_BSEQ_NEXT: {
        return seq_next_item(rhs.at(0), rhs.at(1));
      }
      case yaml::PROD_BSEQ_SCALAR: {
        return any(std::move(rhs.at(3)));
      }
      case yaml::PROD_BSEQ_BSCALAR: {
        return any(std::move(rhs.at(2)));
      }
      case yaml::PROD_BSEQ_BMAP:
      case yaml::PROD_BSEQ_FMAP: {
        return any(std::move(rhs.at(3)));
      }
      case yaml::PROD_BSEQ_BMAP_TRAIL: {
        return any(std::move(rhs.at(4)));
      }
      case yaml::PROD_BSEQ_BSEQ:
      case yaml::PROD_BSEQ_FSEQ: {
        return any(std::move(rhs.at(3)));
      }
      case yaml::PROD_BSEQ_BSEQ_TRAIL: {
        return any(std::move(rhs.at(4)));
      }
      case yaml::PROD_FMAP: {
        return any(std::move(rhs.at(2)));
      }
      case yaml::PROD_FMAP_EMPTY: {
        return InputMap();
      }
      case yaml::PROD_FSEQ: {
        return any(std::move(rhs.at(2)));
      }
      case yaml::PROD_FSEQ_EMPTY: {
        return InputList();
      }
      case yaml::PROD_FSEQ_FIRST: {
        return seq_first_item(rhs.at(0));
      }
      case yaml::PROD_FSEQ_NEXT: {
        return seq_next_item(rhs.at(0), rhs.at(3));
      }
      case yaml::PROD_FSEQ_SCALAR: {
        return any(std::move(rhs.at(1)));
      }
      case yaml::PROD_FSEQ_FSEQ: {
        return any(std::move(rhs.at(1)));
      }
      case yaml::PROD_FSEQ_FMAP: {
        return any(std::move(rhs.at(1)));
      }
      case yaml::PROD_SCALAR_QUOTED:
      case yaml::PROD_MAP_SCALAR_QUOTED: {
        return any(std::move(rhs.at(0)));
      }
      case yaml::PROD_SCALAR_RAW:
      case yaml::PROD_MAP_SCALAR_RAW: {
        std::string text;
        OMEGA_H_CHECK(!rhs.at(0).empty());
        text = any_cast<std::string&&>(std::move(rhs.at(0)));
        text += any_cast<std::string&>(rhs.at(1));
        if (prod == yaml::PROD_MAP_SCALAR_RAW) {
          text += any_cast<std::string&>(rhs.at(2));
        }
        return remove_trailing_whitespace(text);
      }
      case yaml::PROD_SCALAR_HEAD_OTHER:
      case yaml::PROD_SCALAR_HEAD_DOT:
      case yaml::PROD_SCALAR_HEAD_DASH:
      case yaml::PROD_SCALAR_HEAD_DOT_DOT: {
        std::size_t offset;
        if (prod == yaml::PROD_SCALAR_HEAD_OTHER)
          offset = 0;
        else if (prod == yaml::PROD_SCALAR_HEAD_DOT_DOT)
          offset = 2;
        else
          offset = 1;
        char second = any_cast<char>(rhs.at(offset));
        std::string result;
        if (prod == yaml::PROD_SCALAR_HEAD_DOT)
          result += '.';
        else if (prod == yaml::PROD_SCALAR_HEAD_DASH)
          result += '-';
        else if (prod == yaml::PROD_SCALAR_HEAD_DOT_DOT)
          result += "..";
        result += second;
        return any(std::move(result));
      }
      case yaml::PROD_SCALAR_DQUOTED:
      case yaml::PROD_SCALAR_SQUOTED: {
        auto text = any_cast<std::string&&>(std::move(rhs.at(1)));
        text += any_cast<std::string&>(rhs.at(2));
        return any(std::move(text));
      }
      case yaml::PROD_MAP_SCALAR_ESCAPED_EMPTY: {
        return std::string();
      }
      case yaml::PROD_MAP_SCALAR_ESCAPED_NEXT: {
        auto str = any_cast<std::string&&>(std::move(rhs.at(0)));
        str += ',';
        str += any_cast<std::string&>(rhs.at(2));
        return any(std::move(str));
      }
      case yaml::PROD_TAG: {
        return any(std::move(rhs.at(2)));
      }
      case yaml::PROD_BSCALAR: {
        auto parent_indent_level = this->symbol_indentation_stack.at(
            this->symbol_indentation_stack.size() - 5);
        auto& header = any_cast<std::string&>(rhs.at(0));
        auto& leading_empties_or_comments = any_cast<std::string&>(rhs.at(2));
        auto& rest = any_cast<std::string&>(rhs.at(4));
        std::string content;
        std::string ignored_comment;
        handle_block_scalar(parent_indent_level, header,
            leading_empties_or_comments, rest, content, ignored_comment);
        return any(std::move(content));
      }
      case yaml::PROD_BSCALAR_FIRST: {
        return any(std::move(rhs.at(0)));
      }
      // all these cases reduce to concatenating two strings
      case yaml::PROD_BSCALAR_NEXT:
      case yaml::PROD_BSCALAR_LINE:
      case yaml::PROD_DESCAPE_NEXT:
      case yaml::PROD_SESCAPE_NEXT: {
        auto str = any_cast<std::string&&>(std::move(rhs.at(0)));
        str += any_cast<std::string&>(rhs.at(1));
        return any(std::move(str));
      }
      case yaml::PROD_BSCALAR_INDENT: {
        return any(std::move(rhs.at(1)));
      }
      case yaml::PROD_BSCALAR_HEADER_LITERAL:
      case yaml::PROD_BSCALAR_HEADER_FOLDED: {
        std::string result;
        if (prod == yaml::PROD_BSCALAR_HEADER_LITERAL) {
          result += "|";
        } else {
          result += ">";
        }
        auto& rest = any_cast<std::string&>(rhs.at(1));
        result += rest;
        return any(std::move(result));
      }
      case yaml::PROD_DESCAPE: {
        std::string str;
        auto& rest = any_cast<std::string&>(rhs.at(2));
        str += any_cast<char>(rhs.at(1));
        str += rest;
        return any(std::move(str));
      }
      case yaml::PROD_SESCAPE: {
        std::string str;
        auto& rest = any_cast<std::string&>(rhs.at(2));
        str += '\'';
        str += rest;
        return any(std::move(str));
      }
      case yaml::PROD_OTHER_FIRST:
      case yaml::PROD_SPACE_PLUS_FIRST: {
        std::string str;
        str.push_back(any_cast<char>(rhs.at(0)));
        return any(std::move(str));
      }
      case yaml::PROD_SCALAR_TAIL_SPACE:
      case yaml::PROD_SCALAR_TAIL_OTHER:
      case yaml::PROD_DESCAPED_DQUOTED:
      case yaml::PROD_DQUOTED_COMMON:
      case yaml::PROD_SQUOTED_COMMON:
      case yaml::PROD_ANY_COMMON:
      case yaml::PROD_COMMON_SPACE:
      case yaml::PROD_COMMON_OTHER:
      case yaml::PROD_BSCALAR_HEAD_OTHER: {
        return any(std::move(rhs.at(0)));
      }
      // all these cases reduce to appending a character
      case yaml::PROD_DQUOTED_NEXT:
      case yaml::PROD_SQUOTED_NEXT:
      case yaml::PROD_ANY_NEXT:
      case yaml::PROD_SCALAR_TAIL_NEXT:
      case yaml::PROD_SPACE_STAR_NEXT:
      case yaml::PROD_SPACE_PLUS_NEXT:
      case yaml::PROD_BSCALAR_HEAD_NEXT: {
        if (rhs.at(0).empty()) {
          std::stringstream ss;
          ss << "leading characters in " << prod << ": any was empty\n";
          throw ParserFail(ss.str());
        }
        auto str = any_cast<std::string&&>(std::move(rhs.at(0)));
        str += any_cast<char>(rhs.at(1));
        return any(std::move(str));
      }
      case yaml::PROD_DQUOTED_EMPTY:
      case yaml::PROD_SQUOTED_EMPTY:
      case yaml::PROD_ANY_EMPTY:
      case yaml::PROD_DESCAPE_EMPTY:
      case yaml::PROD_SESCAPE_EMPTY:
      case yaml::PROD_SCALAR_TAIL_EMPTY:
      case yaml::PROD_SPACE_STAR_EMPTY:
      case yaml::PROD_BSCALAR_HEAD_EMPTY: {
        return std::string();
      }
      case yaml::PROD_DESCAPED_DQUOT:
      case yaml::PROD_SQUOTED_DQUOT:
      case yaml::PROD_ANY_DQUOT: {
        return '"';
      }
      case yaml::PROD_DESCAPED_SLASH:
      case yaml::PROD_SQUOTED_SLASH:
      case yaml::PROD_ANY_SLASH: {
        return '\\';
      }
      case yaml::PROD_SCALAR_TAIL_SQUOT:
      case yaml::PROD_DQUOTED_SQUOT:
      case yaml::PROD_ANY_SQUOT: {
        return '\'';
      }
      case yaml::PROD_COMMON_COLON: {
        return ':';
      }
      case yaml::PROD_SCALAR_TAIL_DOT:
      case yaml::PROD_COMMON_DOT: {
        return '.';
      }
      case yaml::PROD_SCALAR_TAIL_DASH:
      case yaml::PROD_COMMON_DASH:
      case yaml::PROD_BSCALAR_HEAD_DASH: {
        return '-';
      }
      case yaml::PROD_COMMON_PIPE: {
        return '|';
      }
      case yaml::PROD_COMMON_LSQUARE: {
        return '[';
      }
      case yaml::PROD_COMMON_RSQUARE: {
        return ']';
      }
      case yaml::PROD_COMMON_LCURLY: {
        return '{';
      }
      case yaml::PROD_COMMON_RCURLY: {
        return '}';
      }
      case yaml::PROD_COMMON_RANGLE: {
        return '>';
      }
      case yaml::PROD_COMMON_COMMA: {
        return ',';
      }
      case yaml::PROD_COMMON_PERCENT: {
        return '%';
      }
      case yaml::PROD_COMMON_EXCL: {
        return '!';
      }
    }
    return any();
  }
  any map_first_item(any& first_item) {
    InputMap map;
    OMEGA_H_CHECK(!first_item.empty());
    any map_any = std::move(map);
    return map_next_item(map_any, first_item);
  }
  any map_next_item(any& items, any& next_item) {
    InputMap map = any_cast<InputMap&&>(std::move(items));
    NameValue& pair = any_cast<NameValue&>(next_item);
    map.add(pair.name, std::move(pair.value));
    return any(std::move(map));
  }
  any map_item(any& key_any, any& value_any) {
    NameValue result;
    result.name = any_cast<std::string&&>(std::move(key_any));
    if (value_any.type() == typeid(std::string)) {
      std::string value = any_cast<std::string&&>(std::move(value_any));
      result.value.reset(new InputScalar(value));
    } else if (value_any.type() == typeid(InputList)) {
      InputList value = any_cast<InputList&&>(std::move(value_any));
      result.value.reset(new InputList(std::move(value)));
    } else if (value_any.type() == typeid(InputMap)) {
      InputMap value = any_cast<InputMap&&>(std::move(value_any));
      result.value.reset(new InputMap(std::move(value)));
    } else {
      std::string msg = "unexpected YAML map value type ";
      msg += value_any.type().name();
      msg += " for name \"";
      msg += result.name;
      msg += "\"\n";
      throw ParserFail(msg);
    }
    return any(std::move(result));
  }
  any seq_first_item(any& first_any) {
    InputList list;
    any list_any = std::move(list);
    return seq_next_item(list_any, first_any);
  }
  any seq_next_item(any& items, any& next_item) {
    auto list = any_cast<InputList&&>(std::move(items));
    if (next_item.type() == typeid(std::string)) {
      std::string value = any_cast<std::string&&>(std::move(next_item));
      std::shared_ptr<Input> sptr(new InputScalar(std::move(value)));
      list.add(std::move(sptr));
    } else if (next_item.type() == typeid(InputList)) {
      InputList value = any_cast<InputList&&>(std::move(next_item));
      std::shared_ptr<Input> sptr(new InputList(std::move(value)));
      list.add(std::move(sptr));
    } else if (next_item.type() == typeid(InputMap)) {
      InputMap value = any_cast<InputMap&&>(std::move(next_item));
      std::shared_ptr<Input> sptr(new InputMap(std::move(value)));
      list.add(std::move(sptr));
    } else {
      throw ParserFail(
          "bug in InputYamlReader: unexpected type for sequence item");
    }
    return any(std::move(list));
  }
  /* block scalars are a super complicated mess, this function handles that mess
   */
  void handle_block_scalar(std::size_t parent_indent_level,
      std::string const& header, std::string const& leading_empties_or_comments,
      std::string const& rest, std::string& content, std::string& comment) {
    /* read the header, resulting in: block style, chomping indicator, and
     * indentation indicator */
    char style;
    char chomping_indicator;
    std::size_t indentation_indicator = 0;
    style = header[0];
    std::stringstream ss(header.substr(1, std::string::npos));
    if (header.size() > 1 && my_isdigit(header[1])) {
      ss >> indentation_indicator;
      // indentation indicator is given as a relative number, but we need it in
      // absolute terms
      indentation_indicator += parent_indent_level;
    }
    if (!(ss >> chomping_indicator)) chomping_indicator = '\0';
    /* get information about newlines, indentation level, and comment from
       the leading_empties_or_comments string */
    std::size_t first_newline =
        leading_empties_or_comments.find_first_of("\r\n");
    std::string newline;
    if (first_newline > 0 &&
        leading_empties_or_comments[first_newline - 1] == '\r') {
      newline = "\r\n";
    } else {
      newline = "\n";
    }
    std::size_t keep_beg = first_newline + 1 - newline.size();
    if (leading_empties_or_comments[0] == '#') {
      comment = leading_empties_or_comments.substr(1, keep_beg);
    }
    // according to the YAML spec, a tab is content, not indentation
    std::size_t content_beg =
        leading_empties_or_comments.find_first_not_of("\r\n ");
    if (content_beg == std::string::npos)
      content_beg = leading_empties_or_comments.size();
    std::size_t newline_before_content =
        leading_empties_or_comments.rfind("\n", content_beg);
    std::size_t num_indent_spaces = (content_beg - newline_before_content) - 1;
    /* indentation indicator overrides the derived level of indentation, in case
       the
       user wants to keep some of that indentation as content */
    if (indentation_indicator > 0) {
      if (num_indent_spaces < indentation_indicator) {
        std::string msg = "Indentation indicator ";
        msg += std::to_string(indentation_indicator);
        msg += " > leading spaces ";
        msg += std::to_string(num_indent_spaces);
        msg += "\n";
        throw ParserFail(msg);
      }
      num_indent_spaces = indentation_indicator;
    }
    /* prepend the content from the leading_empties_or_comments to the rest */
    content = leading_empties_or_comments.substr(keep_beg, std::string::npos);
    content += rest;
    /* per Trilinos issue #2090, there can be trailing comments after the block
       scalar which are less indented than it, but they will be included in the
       final NEWLINE token.
       this code removes all contiguous trailing lines which are less indented
       than the content.
     */
    while (true) {
      auto last_newline = content.find_last_of("\n", content.size() - 2);
      if (last_newline == std::string::npos) break;
      std::size_t num_spaces = 0;
      for (auto ispace = last_newline + 1;
           ispace < content.size() && content[ispace] == ' '; ++ispace) {
        ++num_spaces;
      }
      if (num_spaces >= num_indent_spaces) break;
      content.erase(content.begin() + long(last_newline + 1), content.end());
    }
    /* remove both indentation and newlines as dictated by header information */
    std::size_t unindent_pos = 0;
    while (true) {
      std::size_t next_newline = content.find_first_of("\n", unindent_pos);
      if (next_newline == std::string::npos) break;
      std::size_t start_cut = next_newline + 1;
      /* folding block scalars remove newlines */
      if (style == '>') start_cut -= newline.size();
      std::size_t end_cut = next_newline + 1;
      /* the actual amount of indentation in the content varies, start by
         marking it all for removal */
      while (end_cut < content.size() && content[end_cut] == ' ') {
        ++end_cut;
      }
      /* but don't remove more than the actual indent number */
      end_cut = std::min(next_newline + 1 + num_indent_spaces, end_cut);
      /* cut this (newline?)+indentation out of the content */
      content = content.substr(0, start_cut) +
                content.substr(end_cut, std::string::npos);
      unindent_pos = start_cut;
    }
    if (chomping_indicator != '+') {
      content = remove_trailing_whitespace_and_newlines(content);
      if (chomping_indicator != '-') content += newline;
    }
    if (style == '|') {
      // if not already, remove the leading newline
      content = content.substr(newline.size(), std::string::npos);
    }
  }
};

InputYamlReader::~InputYamlReader() {}

static InputMap read_input_without_includes(
    Omega_h::filesystem::path const& path) {
  std::ifstream stream(path.c_str());
  if (!stream.is_open()) {
    Omega_h_fail("Couldn't open Input file \"%s\"\n", path.c_str());
  }
  Omega_h::InputYamlReader reader;
  auto result_any = reader.read_stream(stream, path.string());
  return any_cast<InputMap&&>(std::move(result_any));
}

static bool handle_one_include(InputList& list);

static bool handle_one_include(InputMap& map) {
  for (auto& key : map) {
    if (key == "include") {
      auto const path = map.get<std::string>(key);
      map.remove(key);
      auto included_map = read_input_without_includes(path);
      while (!included_map.map.empty()) {
        auto const it = included_map.map.begin();
        map.add(it->first, std::move(it->second));
        included_map.map.erase(it);
      }
      return true;
    } else if (map.is_list(key)) {
      if (handle_one_include(map.get_list(key))) return true;
    } else if (map.is_map(key)) {
      if (handle_one_include(map.get_map(key))) return true;
    }
  }
  return false;
}

static bool handle_one_include(InputList& list) {
  for (LO i = 0; i < list.size(); ++i) {
    if (list.is_list(i)) {
      if (handle_one_include(list.get_list(i))) return true;
    } else if (list.is_map(i)) {
      if (handle_one_include(list.get_map(i))) return true;
    }
  }
  return false;
}

InputMap read_input(Omega_h::filesystem::path const& path) {
  OMEGA_H_TIME_FUNCTION;
  InputMap map = read_input_without_includes(path);
  while (handle_one_include(map))
    ;
  return map;
}

void update_class_sets(ClassSets* p_sets, InputMap& pl) {
  ClassSets& sets = *p_sets;
  for (auto& set_name : pl) {
    auto& pairs = pl.get_list(set_name);
    auto npairs = pairs.size();
    for (decltype(npairs) i = 0; i < npairs; ++i) {
      auto& pair = pairs.get_list(i);
      if (pair.size() != 2) {
        Omega_h_fail(
            "Expected \"%s\" to be an array of int pairs\n", set_name.c_str());
      }
      auto class_dim = Int(pair.get<int>(0));
      auto class_id = LO(pair.get<int>(1));
      sets[set_name].push_back({class_dim, class_id});
    }
  }
}

static void echo_input_recursive(
    std::ostream& stream, Input& input, Int indent) {
  if (is_type<InputScalar>(input)) {
    auto& scalar = as_type<InputScalar>(input);
    if (std::find(scalar.str.begin(), scalar.str.end(), '\n') !=
        scalar.str.end()) {
      stream << " |\n";
      for (auto c : scalar.str) {
        if (c == '\n') {
          for (Int i = 0; i < indent; ++i) {
            stream << "  ";
          }
        }
        stream << c;
      }
    } else {
      stream << " \'" << scalar.str << '\'';
    }
  } else if (is_type<InputMap>(input)) {
    auto& map = as_type<InputMap>(input);
    stream << '\n';
    for (auto& pair : map.map) {
      auto& name = pair.first;
      for (Int i = 0; i < indent; ++i) {
        stream << "  ";
      }
      stream << name << ":";
      echo_input_recursive(stream, *(pair.second), indent + 1);
    }
  } else if (is_type<InputList>(input)) {
    auto& list = as_type<InputList>(input);
    stream << '\n';
    for (LO i = 0; i < list.size(); ++i) {
      for (Int j = 0; j < indent; ++j) {
        stream << "  ";
      }
      stream << "-";
      echo_input_recursive(stream, *(list.entries[std::size_t(i)]), indent + 1);
    }
  }
}

void echo_input(std::ostream& stream, Input& input) {
  echo_input_recursive(stream, input, 0);
}

void check_unused(Input& input) {
  if (!input.used) {
    auto const full_name = get_full_name(input);
    Omega_h_fail("input \"%s\" was not used!\n", full_name.c_str());
  }
  if (is_type<InputMap>(input)) {
    auto& map = as_type<InputMap>(input);
    for (auto& pair : map.map) {
      check_unused(*(pair.second));
    }
  } else if (is_type<InputList>(input)) {
    auto& list = as_type<InputList>(input);
    for (LO i = 0; i < list.size(); ++i) {
      check_unused(*(list.entries[std::size_t(i)]));
    }
  }
}

#define OMEGA_H_EXPL_INST(InputType)                                           \
  template bool is_type<InputType>(Input&);                                    \
  template InputType& as_type<InputType>(Input&);                              \
  template bool InputMap::is_input<InputType>(std::string const& name);        \
  template InputType& InputMap::use_input<InputType>(std::string const& name); \
  template bool InputList::is_input<InputType>(LO i);                          \
  template InputType& InputList::use_input<InputType>(LO i);
OMEGA_H_EXPL_INST(InputScalar)
OMEGA_H_EXPL_INST(InputMap)
OMEGA_H_EXPL_INST(InputList)
#undef OMEGA_H_EXPL_INST

#define OMEGA_H_EXPL_INST(ScalarType)                                          \
  template ScalarType InputScalar::get<ScalarType>() const;                    \
  template bool InputMap::is<ScalarType>(std::string const& name);             \
  template ScalarType InputMap::get<ScalarType>(std::string const& name);      \
  template ScalarType InputMap::get<ScalarType>(                               \
      std::string const& name, char const* default_value);                     \
  template bool InputList::is<ScalarType>(LO i);                               \
  template ScalarType InputList::get<ScalarType>(LO i);
OMEGA_H_EXPL_INST(std::string)
OMEGA_H_EXPL_INST(bool)
OMEGA_H_EXPL_INST(double)
OMEGA_H_EXPL_INST(int)
OMEGA_H_EXPL_INST(long long)
#undef OMEGA_H_EXPL_INST
}  // namespace Omega_h
