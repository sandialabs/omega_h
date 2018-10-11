
namespace Omega_h {

class YamlInputReader : public Reader {
 public:
  Reader():Reader(yaml::ask_reader_tables()) {}
  ~Reader() override = default;
 protected:
  enum {
    TRIM_NORMAL,
    TRIM_DASH
  };
  virtual void at_shift(any& result_any, int token, std::string& text) {
    using std::swap;
    switch (token) {
      case yaml::TOK_NEWLINE: {
        std::string& result = make_any_ref<std::string>(result_any);
        swap(result, text);
        break;
      }
      case yaml::TOK_SPACE:
      case yaml::TOK_OTHER: {
        result_any = text.at(0);
        break;
      }
    }
  }
  virtual void at_reduce(any& result_any, int prod, std::vector<any>& rhs) {
    using std::swap;
    switch (prod) {
      case yaml::PROD_DOC:
      case yaml::PROD_DOC2: {
        std::size_t offset = prod == yaml::PROD_DOC2 ? 1 : 0;
        TEUCHOS_ASSERT(!rhs.at(offset).empty());
        swap(result_any, rhs.at(offset));
        TEUCHOS_ASSERT(result_any.type() == typeid(ParameterList));
        break;
      }
      case yaml::PROD_TOP_BMAP: {
        TEUCHOS_ASSERT(!rhs.at(0).empty());
        TEUCHOS_ASSERT(rhs.at(0).type() == typeid(PLPair));
        PLPair& pair = any_ref_cast<PLPair>(rhs.at(0));
        any& pair_rhs_any = pair.value.getAny(/* set isUsed = */ false);
        result_any = pair_rhs_any;
        break;
      }
      case yaml::PROD_TOP_FIRST: {
        if (rhs.at(0).type() == typeid(ParameterList)) {
          swap(result_any, rhs.at(0));
        }
        break;
      }
      case yaml::PROD_TOP_NEXT: {
        if (rhs.at(1).type() == typeid(ParameterList)) {
          TEUCHOS_TEST_FOR_EXCEPTION(!rhs.at(0).empty(), ParserFail,
              "Can't specify multiple top-level ParameterLists in one YAML file!\n");
          swap(result_any, rhs.at(1));
        } else {
          swap(result_any, rhs.at(0));
        }
        break;
      }
      case yaml::PROD_BMAP_FIRST:
      case yaml::PROD_FMAP_FIRST: {
        TEUCHOS_ASSERT(rhs.at(0).type() == typeid(PLPair));
        map_first_item(result_any, rhs.at(0));
        TEUCHOS_ASSERT(result_any.type() == typeid(ParameterList));
        break;
      }
      case yaml::PROD_BMAP_NEXT: {
        map_next_item(result_any, rhs.at(0), rhs.at(1));
        break;
      }
      case yaml::PROD_FMAP_NEXT: {
        map_next_item(result_any, rhs.at(0), rhs.at(3));
        break;
      }
      case yaml::PROD_BMAP_SCALAR:
      case yaml::PROD_FMAP_SCALAR:
      case yaml::PROD_FMAP_FMAP:
      case yaml::PROD_FMAP_FSEQ: {
        int scalar_type = interpret_tag(rhs.at(3));
        map_item(result_any, rhs.at(0), rhs.at(4), scalar_type);
        break;
      }
      case yaml::PROD_BMAP_BSCALAR: {
        map_item(result_any, rhs.at(0), rhs.at(3), Scalar::STRING);
        break;
      }
      case yaml::PROD_BMAP_BVALUE: {
        map_item(result_any, rhs.at(0), rhs.at(4));
        break;
      }
      case yaml::PROD_BVALUE_EMPTY: {
        result_any = ParameterList();
        break;
      }
      case yaml::PROD_BVALUE_BMAP:
      case yaml::PROD_BVALUE_BSEQ: {
        swap(result_any, rhs.at(1));
        break;
      }
      case yaml::PROD_BMAP_FMAP: {
        map_item(result_any, rhs.at(0), rhs.at(4));
        break;
      }
      case yaml::PROD_BMAP_FSEQ: {
        TEUCHOS_ASSERT(rhs.at(4).type() == typeid(Array<Scalar>) ||
            rhs.at(4).type() == typeid(Array<Array<Scalar>>));
        int scalar_type = interpret_tag(rhs.at(3));
        map_item(result_any, rhs.at(0), rhs.at(4), scalar_type);
        break;
      }
      case yaml::PROD_BSEQ_FIRST: {
        seq_first_item(result_any, rhs.at(0));
        break;
      }
      case yaml::PROD_BSEQ_NEXT: {
        seq_next_item(result_any, rhs.at(0), rhs.at(1));
        break;
      }
      case yaml::PROD_BSEQ_SCALAR: {
        swap(result_any, rhs.at(3));
        Scalar& scalar = any_ref_cast<Scalar>(result_any);
        scalar.tag_type = interpret_tag(rhs.at(2));
        break;
      }
      case yaml::PROD_BSEQ_BSCALAR: {
        swap(result_any, rhs.at(2));
        break;
      }
      case yaml::PROD_BSEQ_BMAP:
      case yaml::PROD_BSEQ_BMAP_TRAIL:
      case yaml::PROD_BSEQ_FMAP: {
        throw ParserFail("Can't interpret a map inside a sequence as a Teuchos Parameter");
      }
      case yaml::PROD_BSEQ_BSEQ: {
        swap(result_any, rhs.at(3));
        break;
      }
      case yaml::PROD_BSEQ_BSEQ_TRAIL: {
        swap(result_any, rhs.at(4));
        break;
      }
      case yaml::PROD_BSEQ_FSEQ: {
        swap(result_any, rhs.at(3));
        break;
      }
      case yaml::PROD_FMAP: {
        swap(result_any, rhs.at(2));
        break;
      }
      case yaml::PROD_FMAP_EMPTY: {
        result_any = ParameterList();
        break;
      }
      case yaml::PROD_FSEQ: {
        swap(result_any, rhs.at(2));
        TEUCHOS_ASSERT(result_any.type() == typeid(Array<Scalar>) ||
            result_any.type() == typeid(Array<Array<Scalar>>));
        break;
      }
      case yaml::PROD_FSEQ_EMPTY: {
        result_any = Array<Scalar>();
        break;
      }
      case yaml::PROD_FSEQ_FIRST: {
        seq_first_item(result_any, rhs.at(0));
        break;
      }
      case yaml::PROD_FSEQ_NEXT: {
        seq_next_item(result_any, rhs.at(0), rhs.at(3));
        break;
      }
      case yaml::PROD_FSEQ_SCALAR: {
        swap(result_any, rhs.at(1));
        Scalar& scalar = any_ref_cast<Scalar>(result_any);
        scalar.tag_type = interpret_tag(rhs.at(0));
        break;
      }
      case yaml::PROD_FSEQ_FSEQ:
      case yaml::PROD_FSEQ_FMAP: {
        swap(result_any, rhs.at(1));
        break;
      }
      case yaml::PROD_SCALAR_QUOTED:
      case yaml::PROD_MAP_SCALAR_QUOTED: {
        swap(result_any, rhs.at(0));
        break;
      }
      case yaml::PROD_SCALAR_RAW:
      case yaml::PROD_MAP_SCALAR_RAW: {
        Scalar& scalar = make_any_ref<Scalar>(result_any);
        TEUCHOS_ASSERT(!rhs.at(0).empty());
        scalar.text = any_ref_cast<std::string>(rhs.at(0));
        scalar.text += any_ref_cast<std::string>(rhs.at(1));
        if (prod == yaml::PROD_MAP_SCALAR_RAW) {
          scalar.text += any_ref_cast<std::string>(rhs.at(2));
        }
        scalar.text = remove_trailing_whitespace(scalar.text);
        scalar.source = Scalar::RAW;
        scalar.tag_type = -1;
        break;
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
        std::string& result = make_any_ref<std::string>(result_any);
        if (prod == yaml::PROD_SCALAR_HEAD_DOT) result += '.';
        else if (prod == yaml::PROD_SCALAR_HEAD_DASH) result += '-';
        else if (prod == yaml::PROD_SCALAR_HEAD_DOT_DOT) result += "..";
        result += second;
        break;
      }
      case yaml::PROD_SCALAR_DQUOTED:
      case yaml::PROD_SCALAR_SQUOTED: {
        std::string& first = any_ref_cast<std::string>(rhs.at(1));
        std::string& rest = any_ref_cast<std::string>(rhs.at(2));
        Scalar& scalar = make_any_ref<Scalar>(result_any);
        scalar.text += first;
        scalar.text += rest;
        if (prod == yaml::PROD_SCALAR_DQUOTED) {
          scalar.source = Scalar::DQUOTED;
        } else if (prod == yaml::PROD_SCALAR_SQUOTED) {
          scalar.source = Scalar::SQUOTED;
        }
        scalar.tag_type = -1;
        break;
      }
      case yaml::PROD_MAP_SCALAR_ESCAPED_EMPTY: {
        result_any = std::string();
        break;
      }
      case yaml::PROD_MAP_SCALAR_ESCAPED_NEXT: {
        swap(result_any, rhs.at(0));
        std::string& str = any_ref_cast<std::string>(result_any);
        str += ',';
        str += any_ref_cast<std::string>(rhs.at(2));
        break;
      }
      case yaml::PROD_TAG: {
        swap(result_any, rhs.at(2));
        break;
      }
      case yaml::PROD_BSCALAR: {
        std::size_t parent_indent_level =
          this->symbol_indentation_stack.at(
              this->symbol_indentation_stack.size() - 5);
        std::string& header = any_ref_cast<std::string>(rhs.at(0));
        std::string& leading_empties_or_comments =
          any_ref_cast<std::string>(rhs.at(2));
        std::string& rest = any_ref_cast<std::string>(rhs.at(4));
        std::string& content = make_any_ref<std::string>(result_any);
        std::string comment;
        handle_block_scalar(
            parent_indent_level,
            header, leading_empties_or_comments, rest,
            content, comment);
        break;
      }
      case yaml::PROD_BSCALAR_FIRST: {
        swap(result_any, rhs.at(0));
        break;
      }
      // all these cases reduce to concatenating two strings
      case yaml::PROD_BSCALAR_NEXT:
      case yaml::PROD_BSCALAR_LINE:
      case yaml::PROD_DESCAPE_NEXT:
      case yaml::PROD_SESCAPE_NEXT: {
        swap(result_any, rhs.at(0));
        std::string& str = any_ref_cast<std::string>(result_any);
        str += any_ref_cast<std::string>(rhs.at(1));
        break;
      }
      case yaml::PROD_BSCALAR_INDENT: {
        swap(result_any, rhs.at(1));
        break;
      }
      case yaml::PROD_BSCALAR_HEADER_LITERAL:
      case yaml::PROD_BSCALAR_HEADER_FOLDED: {
        std::string& result = make_any_ref<std::string>(result_any);
        if (prod == yaml::PROD_BSCALAR_HEADER_LITERAL) {
          result += "|";
        } else {
          result += ">";
        }
        std::string& rest = any_ref_cast<std::string>(rhs.at(1));
        result += rest;
        break;
      }
      case yaml::PROD_DESCAPE: {
        std::string& str = make_any_ref<std::string>(result_any);
        std::string& rest = any_ref_cast<std::string>(rhs.at(2));
        str += any_cast<char>(rhs.at(1));
        str += rest;
        break;
      }
      case yaml::PROD_SESCAPE: {
        std::string& str = make_any_ref<std::string>(result_any);
        std::string& rest = any_ref_cast<std::string>(rhs.at(2));
        str += '\'';
        str += rest;
        break;
      }
      case yaml::PROD_OTHER_FIRST:
      case yaml::PROD_SPACE_PLUS_FIRST: {
        std::string& str = make_any_ref<std::string>(result_any);
        str.push_back(any_cast<char>(rhs.at(0)));
        break;
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
        swap(result_any, rhs.at(0));
        break;
      }
      // all these cases reduce to appending a character
      case yaml::PROD_DQUOTED_NEXT:
      case yaml::PROD_SQUOTED_NEXT:
      case yaml::PROD_ANY_NEXT:
      case yaml::PROD_SCALAR_TAIL_NEXT:
      case yaml::PROD_SPACE_STAR_NEXT:
      case yaml::PROD_SPACE_PLUS_NEXT:
      case yaml::PROD_BSCALAR_HEAD_NEXT: {
        TEUCHOS_TEST_FOR_EXCEPTION(rhs.at(0).empty(), ParserFail,
            "leading characters in " << prod << ": any was empty\n");
        swap(result_any, rhs.at(0));
        std::string& str = any_ref_cast<std::string>(result_any);
        str += any_cast<char>(rhs.at(1));
        break;
      }
      case yaml::PROD_DQUOTED_EMPTY:
      case yaml::PROD_SQUOTED_EMPTY:
      case yaml::PROD_ANY_EMPTY:
      case yaml::PROD_DESCAPE_EMPTY:
      case yaml::PROD_SESCAPE_EMPTY:
      case yaml::PROD_SCALAR_TAIL_EMPTY:
      case yaml::PROD_SPACE_STAR_EMPTY:
      case yaml::PROD_BSCALAR_HEAD_EMPTY: {
        result_any = std::string();
        break;
      }
      case yaml::PROD_DESCAPED_DQUOT:
      case yaml::PROD_SQUOTED_DQUOT:
      case yaml::PROD_ANY_DQUOT: {
        result_any = '"';
        break;
      }
      case yaml::PROD_DESCAPED_SLASH:
      case yaml::PROD_SQUOTED_SLASH:
      case yaml::PROD_ANY_SLASH: {
        result_any = '\\';
        break;
      }
      case yaml::PROD_SCALAR_TAIL_SQUOT:
      case yaml::PROD_DQUOTED_SQUOT:
      case yaml::PROD_ANY_SQUOT: {
        result_any = '\'';
        break;
      }
      case yaml::PROD_COMMON_COLON: {
        result_any = ':';
        break;
      }
      case yaml::PROD_SCALAR_TAIL_DOT:
      case yaml::PROD_COMMON_DOT: {
        result_any = '.';
        break;
      }
      case yaml::PROD_SCALAR_TAIL_DASH:
      case yaml::PROD_COMMON_DASH:
      case yaml::PROD_BSCALAR_HEAD_DASH: {
        result_any = '-';
        break;
      }
      case yaml::PROD_COMMON_PIPE: {
        result_any = '|';
        break;
      }
      case yaml::PROD_COMMON_LSQUARE: {
        result_any = '[';
        break;
      }
      case yaml::PROD_COMMON_RSQUARE: {
        result_any = ']';
        break;
      }
      case yaml::PROD_COMMON_LCURLY: {
        result_any = '{';
        break;
      }
      case yaml::PROD_COMMON_RCURLY: {
        result_any = '}';
        break;
      }
      case yaml::PROD_COMMON_RANGLE: {
        result_any = '>';
        break;
      }
      case yaml::PROD_COMMON_COMMA: {
        result_any = ',';
        break;
      }
      case yaml::PROD_COMMON_PERCENT: {
        result_any = '%';
        break;
      }
      case yaml::PROD_COMMON_EXCL: {
        result_any = '!';
        break;
      }
    }
  }
  void map_first_item(any& result_any, any& first_item) {
    ParameterList& list = make_any_ref<ParameterList>(result_any);
    TEUCHOS_ASSERT(!first_item.empty());
    PLPair& pair = any_ref_cast<PLPair>(first_item);
    safe_set_entry(list, pair.key, pair.value);
  }
  void map_next_item(any& result_any, any& items, any& next_item) {
    using std::swap;
    swap(result_any, items);
    ParameterList& list = any_ref_cast<ParameterList>(result_any);
    PLPair& pair = any_ref_cast<PLPair>(next_item);
    safe_set_entry(list, pair.key, pair.value);
  }
  void map_item(any& result_any, any& key_any, any& value_any, int scalar_type = -1) {
    using std::swap;
    PLPair& result = make_any_ref<PLPair>(result_any);
    {
      std::string& key = any_ref_cast<Scalar>(key_any).text;
      swap(result.key, key);
    }
    resolve_map_value(value_any, scalar_type);
    if (value_any.type() == typeid(bool)) {
      bool value = any_cast<bool>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(int)) {
      int value = any_cast<int>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(long long)) {
      long long value = any_cast<long long>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(double)) {
      double value = any_cast<double>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(std::string)) {
      std::string& value = any_ref_cast<std::string >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<int>)) {
      Array<int>& value = any_ref_cast<Array<int> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<long long>)) {
      Array<long long>& value = any_ref_cast<Array<long long> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<double>)) {
      Array<double>& value = any_ref_cast<Array<double> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<std::string>)) {
      Array<std::string>& value = any_ref_cast<Array<std::string> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(TwoDArray<int>)) {
      TwoDArray<int>& value = any_ref_cast<TwoDArray<int> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(TwoDArray<long long>)) {
      TwoDArray<long long>& value = any_ref_cast<TwoDArray<long long> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(TwoDArray<double>)) {
      TwoDArray<double>& value = any_ref_cast<TwoDArray<double> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(TwoDArray<std::string>)) {
      TwoDArray<std::string>& value = any_ref_cast<TwoDArray<std::string> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(ParameterList)) {
      ParameterList& value = any_ref_cast<ParameterList>(value_any);
      ParameterList& result_pl = result.value.setList();
      swap(result_pl, value);
      result_pl.setName(result.key);
    } else {
      std::string msg = "unexpected YAML map value type ";
      msg += value_any.type().name();
      msg += " for key \"";
      msg += result.key;
      msg += "\"\n";
      throw ParserFail(msg);
    }
  }
  void resolve_map_value(any& value_any, int scalar_type = -1) const {
    if (value_any.type() == typeid(Scalar)) {
      Scalar& scalar_value = any_ref_cast<Scalar>(value_any);
      if (scalar_type == -1) {
        scalar_type = scalar_value.infer_type();
      }
      if (scalar_type == Scalar::BOOL) {
        value_any = parse_as<bool>(scalar_value.text);
      } else if (scalar_type == Scalar::INT) {
        value_any = parse_as<int>(scalar_value.text);
      } else if (scalar_type == Scalar::LONG_LONG) {
        value_any = parse_as<long long>(scalar_value.text);
      } else if (scalar_type == Scalar::DOUBLE) {
        value_any = parse_as<double>(scalar_value.text);
      } else {
        value_any = scalar_value.text;
      }
    } else if (value_any.type() == typeid(Array<Scalar>)) {
      Array<Scalar>& scalars = any_ref_cast<Array<Scalar> >(value_any);
      if (scalar_type == -1) {
        if (scalars.size() == 0) {
          throw ParserFail("implicitly typed arrays can't be empty\n"
                           "(need to determine their element type)\n");
        }
        /* Teuchos::Array uses std::vector but doesn't account for std::vector<bool>,
           so it can't store bools */
        scalar_type = Scalar::INT;
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          scalar_type = std::min(scalar_type, scalars[i].infer_type());
        }
      }
      if (scalar_type == Scalar::INT) {
        Array<int> result(scalars.size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          result[i] = parse_as<int>(scalars[i].text);
        }
        value_any = result;
      } else if (scalar_type == Scalar::LONG_LONG) {
        Array<long long> result(scalars.size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          result[i] = parse_as<long long>(scalars[i].text);
        }
        value_any = result;
      } else if (scalar_type == Scalar::DOUBLE) {
        Array<double> result(scalars.size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          result[i] = parse_as<double>(scalars[i].text);
        }
        value_any = result;
      } else if (scalar_type == Scalar::STRING) {
        Array<std::string> result(scalars.size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          result[i] = scalars[i].text;
        }
        value_any = result;
      }
    } else if (value_any.type() == typeid(Array<Array<Scalar>>)) {
      Array<Array<Scalar>>& scalars = any_ref_cast<Array<Array<Scalar>> >(value_any);
      if (scalar_type == -1) {
        if (scalars.size() == 0) {
          throw ParserFail("implicitly typed 2D arrays can't be empty\n"
                           "(need to determine their element type)\n");
        }
        /* Teuchos::Array uses std::vector but doesn't account for std::vector<bool>,
           so it can't store bools */
        scalar_type = Scalar::INT;
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          if (scalars[0].size() == 0) {
            throw ParserFail("implicitly typed 2D arrays can't have empty rows\n"
                             "(need to determine their element type)\n");
          }
          if (scalars[i].size() != scalars[0].size()) {
            throw ParserFail("2D array: sub-arrays are different sizes");
          }
          for (Teuchos_Ordinal j = 0; j < scalars[i].size(); ++j) {
            int item_type = scalars[i][j].infer_type();
            scalar_type = std::min(scalar_type, item_type);
          }
        }
      }
      if (scalar_type == Scalar::INT) {
        TwoDArray<int> result(scalars.size(), scalars[0].size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          for (Teuchos_Ordinal j = 0; j < scalars[0].size(); ++j) {
            result(i, j) = parse_as<int>(scalars[i][j].text);
          }
        }
        value_any = result;
      } else if (scalar_type == Scalar::LONG_LONG) {
        TwoDArray<long long> result(scalars.size(), scalars[0].size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          for (Teuchos_Ordinal j = 0; j < scalars[0].size(); ++j) {
            result(i, j) = parse_as<long long>(scalars[i][j].text);
          }
        }
        value_any = result;
      } else if (scalar_type == Scalar::DOUBLE) {
        TwoDArray<double> result(scalars.size(), scalars[0].size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          for (Teuchos_Ordinal j = 0; j < scalars[0].size(); ++j) {
            result(i, j) = parse_as<double>(scalars[i][j].text);
          }
        }
        value_any = result;
      } else if (scalar_type == Scalar::STRING) {
        TwoDArray<std::string> result(scalars.size(), scalars[0].size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          for (Teuchos_Ordinal j = 0; j < scalars[0].size(); ++j) {
            result(i, j) = scalars[i][j].text;
          }
        }
        value_any = result;
      }
    }
  }
  int interpret_tag(any& tag_any) {
    if (tag_any.type() != typeid(std::string)) return -1;
    std::string& text = any_ref_cast<std::string>(tag_any);
    if (text.find("bool") != std::string::npos) return Scalar::BOOL;
    else if (text.find("int") != std::string::npos) return Scalar::INT;
    else if (text.find("double") != std::string::npos) return Scalar::DOUBLE;
    else if (text.find("string") != std::string::npos) return Scalar::STRING;
    else {
      std::string msg = "Unable to parse type tag \"";
      msg += text;
      msg += "\"\n";
      throw ParserFail(msg);
    }
  }
  void seq_first_item(any& result_any, any& first_any) {
    using std::swap;
    if (first_any.type() == typeid(Scalar)) {
      Array<Scalar>& a = make_any_ref<Array<Scalar> >(result_any);
      Scalar& v = any_ref_cast<Scalar>(first_any);
      a.push_back(Scalar());
      swap(a.back(), v);
    } else if (first_any.type() == typeid(Array<Scalar>)) {
      Array<Array<Scalar>>& a = make_any_ref<Array<Array<Scalar>> >(result_any);
      Array<Scalar>& v = any_ref_cast<Array<Scalar> >(first_any);
      a.push_back(Array<Scalar>());
      swap(a.back(), v);
    } else {
      throw Teuchos::ParserFail(
          "bug in YAMLParameterList::Reader: unexpected type for first sequence item");
    }
  }
  void seq_next_item(any& result_any, any& items, any& next_item) {
    using std::swap;
    swap(result_any, items);
    if (result_any.type() == typeid(Array<Scalar>)) {
      Array<Scalar>& a = any_ref_cast<Array<Scalar> >(result_any);
      Scalar& val = any_ref_cast<Scalar>(next_item);
      a.push_back(Scalar());
      swap(a.back(), val);
    } else if (result_any.type() == typeid(Array<Array<Scalar>>)) {
      Array<Array<Scalar>>& a = any_ref_cast<Array<Array<Scalar>> >(result_any);
      Array<Scalar>& v = any_ref_cast<Array<Scalar> >(next_item);
      a.push_back(Array<Scalar>());
      swap(a.back(), v);
    } else {
      throw Teuchos::ParserFail(
          "bug in YAMLParameterList::Reader: unexpected type for next sequence item");
    }
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
      TEUCHOS_TEST_FOR_EXCEPTION(num_indent_spaces < indentation_indicator,
          Teuchos::ParserFail,
          "Indentation indicator " << indentation_indicator << " > leading spaces " << num_indent_spaces);
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
}
