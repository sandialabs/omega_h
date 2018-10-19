#include <Omega_h_input.hpp>

namespace Omega_h {

Input::Input() : parent(nullptr), used(false) {}

std::string get_full_name(Input& input) {
  std::string full_name;
  if (input.parent != nullptr) {
    auto& parent = *(input.parent);
    full_name = get_full_name(parent);
    if (is_type<InputList>(parent)) {
      auto i = as_type<InputList>(parent).position(input);
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

template <class InputType>
bool is_type<InputType>(Input& input) {
  auto ptr_in = &input;
  auto ptr_out = dynamic_cast<InputType*>(ptr_in);
  return ptr_out != nullptr;
}

template <class InputType>
InputType& as_type<InputType>(Input& input) {
  return dynamic_cast<InputType&>(input);
}

InputScalar::InputScalar(std::string const& str_in) : str(str_in) {}

bool InputScalar::as(std::string& out) { out = str; return true; }

bool InputScalar::as(bool& out) {
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

bool InputScalar::as(double& out) {
  std::istringstream ss(str);
  double val;
  ss >> std::noskipws >> val;
  return ss.eof() && !ss.fail();
}

bool InputScalar::as(int& out) {
  std::istringstream ss(text);
  using LL = long long;
  LL val;
  ss >> std::noskipws >> val;
  if (ss.eof() && !ss.fail() &&
    (val >= LL(std::numeric_limits<int>::min())) &&
    (val <= LL(std::numeric_limits<int>::max()))) {
    out = int(val);
    return true;
  }
  return false;
}

bool InputScalar::as(long long& out) {
  std::istringstream ss(str);
  long long val;
  ss >> std::noskipws >> val;
  return ss.eof() && !ss.fail();
}

template <class T>
T InputScalar::get() const {
  T out;
  if (!as(out)) {
    auto full_name = get_full_name(*this);
    Omega_h_fail("InputScalar \"%s\" string \"%s\" is not interpretable as a %s",
        full_name.c_str(),
        str.c_str(),
        (std::is_same<T, int>::value ? "int" :
        (std::is_same<T, double>::value ? "double" :
        (std::is_same<T, long long>::value ? "long long" : "unknown type")))
        );
  }
  return out;
}


void InputMap::add(std::string const& name,
    std::unique_ptr<Input>&& input) {
  auto it = map.upper_bound(name);
  input->parent = this;
  auto const did = entries.emplace(name, std::move(input)).second;
  if (!did) {
    fail(
        "tried to add already existing InputMap name \"%s\"\n",
        name.c_str());
  }
}

template <class InputType>
bool InputMap::is_input(std::string const& name) {
  auto const it = map.find(name);
  if (it == map.end()) return false;
  auto const& uptr = *it;
  return is_type<InputType>(*uptr);
}

template <class ScalarType>
bool InputMap::is(std::string const& name) {
  auto const it = map.find(name);
  if (it == map.end()) return false;
  auto const& uptr = *it;
  ScalarType ignored;
  return is_type<InputScalar>(*uptr) &&
    as_type<InputScalar>(*uptr).as(ignored);
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
    auto s = get_full_name(*this) + name;
    fail("tried to find InputMap entry \"%s\" that doesn't exist\n", s.c_str());
  }
  return *(it->second);
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
  return this->use_input<InputMap>(name);
}

InputList& InputMap::get_list(std::string const& name) {
  return this->use_input<InputList>(name);
}

template <class ScalarType>
ScalarType& InputMap::get(std::string const& name, char const* default_value) {
  if (this->is<ScalarType>(name)) {
    return this->get<ScalarType>(name);
  }
  this->add(name, InputScalar(default_value));
}

void InputList::add(std::unique_ptr<Input>&& input) {
  input->parent = this;
  entries.push_back(std::move(input));
}

LO InputList::position(Input& input) {
  auto it = std::find(entries.begin(), entries.end(),
      [&](std::unique_ptr<Input> const& uptr) {
        return uptr.get() == &input;
      });
  OMEGA_H_CHECK(it != entries.end());
  return LO(it - entries.begin());
}

LO InputList::size() { return LO(entries.size()); }

Input& InputList::at(LO i) {
  return *(entries.at(std::size_t(i)));
}

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
  auto const& input = at(i);
  ScalarType ignored;
  return is_type<InputScalar>(input) &&
    as_type<InputScalar>(input).as(ignored);
}

bool InputList::is_map(LO i) {
  return this->is_input<InputMap>(i);
}

bool InputList::is_list(LO i) {
  return this->is_input<InputList>(i);
}

template <class ScalarType>
ScalarType InputList::get(LO i) {
  return this->use_input<InputScalar>(i).get<ScalarType>();
}

InputMap& InputList::get_map(LO i) {
  return this->use_input<InputMap>(i);
}

InputList& InputList::get_list(LO i) {
  return this->use_input<InputList>(i);
}

struct NameValue {
  std::string name;
  std::unique_ptr<Input> value;
};

class InputYamlReader : public Reader {
 public:
  Reader():Reader(yaml::ask_reader_tables()) {}
  ~Reader() override = default;
 protected:
  enum {
    TRIM_NORMAL,
    TRIM_DASH
  };
  virtual any at_shift(int token, std::string& text) {
    switch (token) {
      case yaml::TOK_NEWLINE: {
        return std::move(text);
      }
      case yaml::TOK_SPACE:
      case yaml::TOK_OTHER: {
        return text.at(0);
      }
    }
  }
  virtual any at_reduce(int prod, std::vector<any>& rhs) {
    switch (prod) {
      case yaml::PROD_DOC:
      case yaml::PROD_DOC2: {
        std::size_t offset = prod == yaml::PROD_DOC2 ? 1 : 0;
        OMEGA_H_CHECK(!rhs.at(offset).empty());
        auto& result_any = rhs.at(offset);
        OMEGA_H_CHECK(result_any.type() == typeid(InputMap));
        return std::move(result_any);
      }
      case yaml::PROD_TOP_BMAP: {
        OMEGA_H_CHECK(!rhs.at(0).empty());
        OMEGA_H_CHECK(rhs.at(0).type() == typeid(NameValue));
        return std::move(as_type<InputMap>(*(any_cast<NameValue&>(rhs.at(0)).value)));
      }
      case yaml::PROD_TOP_FIRST: {
        if (rhs.at(0).type() == typeid(InputMap)) {
          return std::move(rhs.at(0));
        }
        return any();
      }
      case yaml::PROD_TOP_NEXT: {
        if (rhs.at(1).type() == typeid(InputMap)) {
          if (!rhs.at(0).empty()) {
            throw ParserFail("Can't specify multiple top-level InputMaps in one YAML file!\n");
          }
          return std::move(rhs.at(1));
        } else {
          return std::move(rhs.at(0));
        }
      }
      case yaml::PROD_BMAP_FIRST:
      case yaml::PROD_FMAP_FIRST: {
        OMEGA_H_CHECK(rhs.at(0).type() == typeid(NameValue));
        auto result_any = map_first_item(result_any, rhs.at(0));
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
      case yaml::PROD_FMAP_SCALAR {
        auto& scalar = any_cast<InputScalar&>(rhs.at(4));
        std::unique_ptr<Input> uptr(new InputScalar(std::move(scalar));
        any item = std::move(uptr);
        return map_item(rhs.at(0), item);
      }
      case yaml::PROD_FMAP_FMAP: {
        auto& map = any_cast<InputMap&>(rhs.at(4));
        std::unique_ptr<Input> uptr(new InputMap(std::move(map));
        any item = std::move(uptr);
        return map_item(rhs.at(0), item);
      }
      case yaml::PROD_FMAP_FSEQ: {
        auto& list = any_cast<InputMap&>(rhs.at(4));
        std::unique_ptr<Input> uptr(new InputList(std::move(list));
        any item = std::move(uptr);
        return map_item(rhs.at(0), item);
      }
      case yaml::PROD_BMAP_BSCALAR: {
        return map_item(rhs.at(0), rhs.at(3));
      }
      case yaml::PROD_BMAP_BVALUE: {
        return map_item(rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_BVALUE_EMPTY: {
        std::unique_ptr<Input> uptr(new InputMap());
        return uptr;
      }
      case yaml::PROD_BVALUE_BMAP: {
        auto& map = any_cast<InputMap&>(rhs.at(1));
        std::unique_ptr<Input> uptr(new InputMap(std::move(map));
        return uptr;
      }
      case yaml::PROD_BVALUE_BSEQ: {
        auto& list = any_cast<InputList&>(rhs.at(1));
        std::unique_ptr<Input> uptr(new InputList(std::move(list));
        return uptr;
      }
      case yaml::PROD_BMAP_FMAP: {
        return map_item(rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_BMAP_FSEQ: {
        return map_item(result_any, rhs.at(0), rhs.at(4));
      }
      case yaml::PROD_BSEQ_FIRST: {
        return seq_first_item(rhs.at(0));
      }
      case yaml::PROD_BSEQ_NEXT: {
        return seq_next_item(rhs.at(0), rhs.at(1));
      }
      case yaml::PROD_BSEQ_SCALAR: {
        auto& scalar = any_cast<InputScalar&>(rhs.at(4));
        std::unique_ptr<Input> uptr(new InputScalar(std::move(scalar));
        return uptr;
      }
      case yaml::PROD_BSEQ_BSCALAR: {
        auto& scalar = any_cast<InputScalar&>(rhs.at(2));
        std::unique_ptr<Input> uptr(new InputScalar(std::move(scalar));
        return uptr;
      }
      case yaml::PROD_BSEQ_BMAP:
      case yaml::PROD_BSEQ_FMAP: {
        auto& map = any_cast<InputMap&>(rhs.at(3));
        std::unique_ptr<Input> uptr(new InputMap(std::move(map));
        return uptr;
      }
      case yaml::PROD_BSEQ_BMAP_TRAIL: {
        auto& map = any_cast<InputMap&>(rhs.at(4));
        std::unique_ptr<Input> uptr(new InputMap(std::move(map));
        return uptr;
      }
      case yaml::PROD_BSEQ_BSEQ:
      case yaml::PROD_BSEQ_FSEQ: {
        auto& list = any_cast<InputList&>(rhs.at(3));
        std::unique_ptr<Input> uptr(new InputList(std::move(list));
        return uptr;
      }
      case yaml::PROD_BSEQ_BSEQ_TRAIL: {
        auto& list = any_cast<InputList&>(rhs.at(4));
        std::unique_ptr<Input> uptr(new InputList(std::move(list));
        return uptr;
      }
      case yaml::PROD_FMAP: {
        return std::move(rhs.at(2));
      }
      case yaml::PROD_FMAP_EMPTY: {
        return InputMap();
      }
      case yaml::PROD_FSEQ: {
        return std::move(rhs.at(2));
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
        auto& scalar = any_cast<InputScalar&>(rhs.at(1));
        std::unique_ptr<Input> uptr(new InputScalar(std::move(scalar));
        return uptr;
      }
      case yaml::PROD_FSEQ_FSEQ: {
        auto& list = any_cast<InputList&>(rhs.at(1));
        std::unique_ptr<Input> uptr(new InputList(std::move(list));
        return uptr;
      }
      case yaml::PROD_FSEQ_FMAP: {
        auto& map = any_cast<InputMap&>(rhs.at(1));
        std::unique_ptr<Input> uptr(new InputMap(std::move(map));
        return uptr;
      }
      case yaml::PROD_SCALAR_QUOTED:
      case yaml::PROD_MAP_SCALAR_QUOTED: {
        return std::move(rhs.at(0));
      }
      case yaml::PROD_SCALAR_RAW:
      case yaml::PROD_MAP_SCALAR_RAW: {
        std::string text;
        OMEGA_H_CHECK(!rhs.at(0).empty());
        text = any_cast<std::string&&>(rhs.at(0));
        text += any_cast<std::string&>(rhs.at(1));
        if (prod == yaml::PROD_MAP_SCALAR_RAW) {
          text += any_cast<std::string&>(rhs.at(2));
        }
        text = remove_trailing_whitespace(scalar.str);
        return text;
      }
      case yaml::PROD_SCALAR_HEAD_OTHER:
      case yaml::PROD_SCALAR_HEAD_DOT:
      case yaml::PROD_SCALAR_HEAD_DASH:
      case yaml::PROD_SCALAR_HEAD_DOT_DOT: {
        std::size_t offset;
        if (prod == yaml::PROD_SCALAR_HEAD_OTHER) offset = 0;
        else if (prod == yaml::PROD_SCALAR_HEAD_DOT_DOT) offset = 2;
        else offset = 1;
        char second = any_cast<char>(rhs.at(offset));
        std::string result;
        if (prod == yaml::PROD_SCALAR_HEAD_DOT) result += '.';
        else if (prod == yaml::PROD_SCALAR_HEAD_DASH) result += '-';
        else if (prod == yaml::PROD_SCALAR_HEAD_DOT_DOT) result += "..";
        result += second;
        return result;
      }
      case yaml::PROD_SCALAR_DQUOTED:
      case yaml::PROD_SCALAR_SQUOTED: {
        auto text = any_cast<std::string&&>(rhs.at(1));
        text += any_cast<std::string&>(rhs.at(2));
        return text;
      }
      case yaml::PROD_MAP_SCALAR_ESCAPED_EMPTY: {
        return std::string();
      }
      case yaml::PROD_MAP_SCALAR_ESCAPED_NEXT: {
        auto str = any_cast<std::string&&>(rhs.at(0));
        str += ',';
        str += any_cast<std::string&>(rhs.at(2));
        return std::move(str);
      }
      case yaml::PROD_TAG: {
        return std::move(rhs.at(2));
      }
      case yaml::PROD_BSCALAR: {
        auto parent_indent_level =
          this->symbol_indentation_stack.at(
              this->symbol_indentation_stack.size() - 5);
        auto& header = any_cast<std::string&>(rhs.at(0));
        auto& leading_empties_or_comments =
          any_cast<std::string&>(rhs.at(2));
        auto& rest = any_cast<std::string&>(rhs.at(4));
        std::string content;
        std::string ignored_comment;
        handle_block_scalar(
            parent_indent_level,
            header, leading_empties_or_comments, rest,
            content, ignored_comment);
        return content;
      }
      case yaml::PROD_BSCALAR_FIRST: {
        return std::move(rhs.at(0));
      }
      // all these cases reduce to concatenating two strings
      case yaml::PROD_BSCALAR_NEXT:
      case yaml::PROD_BSCALAR_LINE:
      case yaml::PROD_DESCAPE_NEXT:
      case yaml::PROD_SESCAPE_NEXT: {
        auto str = any_cast<std::string&&>(rhs.at(0));
        str += any_cast<std::string&>(rhs.at(1));
        return str;
      }
      case yaml::PROD_BSCALAR_INDENT: {
        return std::move(rhs.at(1));
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
        return result;
      }
      case yaml::PROD_DESCAPE: {
        std::string str;
        auto& rest = any_cast<std::string&>(rhs.at(2));
        str += any_cast<char>(rhs.at(1));
        str += rest;
        return str;
      }
      case yaml::PROD_SESCAPE: {
        std::string str;
        auto& rest = any_cast<std::string&>(rhs.at(2));
        str += '\'';
        str += rest;
        return str;
      }
      case yaml::PROD_OTHER_FIRST:
      case yaml::PROD_SPACE_PLUS_FIRST: {
        std::string str;
        str.push_back(any_cast<char>(rhs.at(0)));
        return str;
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
        return std::move(rhs.at(0));
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
          ss << "leading characters in " << prod << ": any was empty\n");
          throw ParserFail(ss.str());
        }
        auto str = any_cast<std::string&&>(rhs.at(0));
        str += any_cast<char>(rhs.at(1));
        return str;
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
  }
  any map_first_item(any& first_item) {
    InputMap map;
    OMEGA_H_CHECK(!first_item.empty());
    any map_any = map;
    return map_next_item(map_any, first_item);
  }
  any map_next_item(any& items, any& next_item) {
    InputMap map = any_cast<InputMap&&>(items);
    NameValue& pair = any_cast<NameValue&>(next_item);
    map.add(pair.name, std::move(pair.value));
    return map;
  }
  any map_item(any& key_any, any& value_any) {
    NameValue result;
    result.key = any_cast<std::string&&>(key_any);
    if (value_any.type() == typeid(std::string)) {
      std::string value = any_cast<std::string&&>(value_any);
      result.value.reset(new InputScalar(value));
    } else if (value_any.type() == typeid(InputList)) {
      InputList value = any_cast<InputList&&>(value_any);
      result.value.reset(new InputList(std::move(value)));
    } else if (value_any.type() == typeid(ParameterList)) {
      InputMap value = any_cast<InputMap&&>(value_any);
      result.value.reset(new InputMap(std::move(value)));
    } else {
      std::string msg = "unexpected YAML map value type ";
      msg += value_any.type().name();
      msg += " for key \"";
      msg += result.key;
      msg += "\"\n";
      throw ParserFail(msg);
    }
    return result;
  }
  any seq_first_item(any& first_any) {
    InputList list;
    any list_any = list;
    return seq_next_item(list, first_any);
  }
  any seq_next_item(any& items, any& next_item) {
    auto list = any_cast<InputList&&>(items);
    if (first_any.type() == typeid(std::string)) {
      std::string value = any_cast<std::string&&>(value_any);
      list.add(new InputScalar(std::move(value)));
    } else if (first_any.type() == typeid(InputList)) {
      InputList value = any_cast<InputList&&>(value_any);
      list.add(new InputList(std::move(value)));
    } else if (first_any.type() == typeid(InputMap)) {
      InputMap value = any_cast<InputMap&&>(value_any);
      list.add(new InputMap(std::move(value)));
    } else {
      throw ParserFail(
          "bug in YAMLParameterList::Reader: unexpected type for sequence item");
    }
    return list;
  }
  /* block scalars are a super complicated mess, this function handles that mess */
  void handle_block_scalar(
      std::size_t parent_indent_level,
      std::string const& header,
      std::string const& leading_empties_or_comments,
      std::string const& rest,
      std::string& content,
      std::string& comment) {
    /* read the header, resulting in: block style, chomping indicator, and indentation indicator */
    char style;
    char chomping_indicator;
    std::size_t indentation_indicator = 0;
    style = header[0];
    std::stringstream ss(header.substr(1,std::string::npos));
    if (header.size() > 1 && my_isdigit(header[1])) {
      ss >> indentation_indicator;
      // indentation indicator is given as a relative number, but we need it in absolute terms
      indentation_indicator += parent_indent_level;
    }
    if (!(ss >> chomping_indicator)) chomping_indicator = '\0';
    /* get information about newlines, indentation level, and comment from
       the leading_empties_or_comments string */
    std::size_t first_newline = leading_empties_or_comments.find_first_of("\r\n");
    std::string newline;
    if (first_newline > 0 && leading_empties_or_comments[first_newline - 1] == '\r') {
      newline = "\r\n";
    } else {
      newline = "\n";
    }
    std::size_t keep_beg = first_newline + 1 - newline.size();
    if (leading_empties_or_comments[0] == '#') {
      comment = leading_empties_or_comments.substr(1, keep_beg);
    }
    // according to the YAML spec, a tab is content, not indentation
    std::size_t content_beg = leading_empties_or_comments.find_first_not_of("\r\n ");
    if (content_beg == std::string::npos) content_beg = leading_empties_or_comments.size();
    std::size_t newline_before_content = leading_empties_or_comments.rfind("\n", content_beg);
    std::size_t num_indent_spaces = (content_beg - newline_before_content) - 1;
    /* indentation indicator overrides the derived level of indentation, in case the
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
           ispace < content.size() && content[ispace] == ' ';
           ++ispace) {
        ++num_spaces;
      }
      if (num_spaces >= num_indent_spaces) break;
      content.erase(content.begin() + last_newline + 1, content.end());
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

#define OMEGA_H_EXPL_INST(InputType) \
template bool is_type<InputType>(Input&); \
template InputType as_type<InputType>(Input&); \
template bool InputMap::is_input<InputType>(std::string const& name); \
template InputType& InputMap::use_input<InputType>(std::string const& name) \
template bool InputList::is_input(LO i); \
template InputType& InputList::use_input(LO i);
OMEGA_H_EXPL_INST_INPUT(InputScalar)
OMEGA_H_EXPL_INST_INPUT(InputMap)
OMEGA_H_EXPL_INST_INPUT(InputList)
#undef OMEGA_H_EXPL_INST

#define OMEGA_H_EXPL_INST(ScalarType) \
template ScalarType InputScalar::get<ScalarType>() const; \
template bool InputMap::is<ScalarType>(std::string const& name); \
template ScalarType InputMap::get<ScalarType>(std::string const& name); \
template ScalarType& InputMap::get(std::string const& name, char const* default_value); \
template bool InputList::is<ScalarType>(LO i); \
template ScalarType InputList::get<ScalarType>(LO i);
OMEGA_H_EXPL_INST(std::string)
OMEGA_H_EXPL_INST(double)
OMEGA_H_EXPL_INST(int)
OMEGA_H_EXPL_INST(long long)
#undef OMEGA_H_EXPL_INST

}
