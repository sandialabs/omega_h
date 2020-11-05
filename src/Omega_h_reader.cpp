#include "Omega_h_reader.hpp"

#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>

#include "Omega_h_fail.hpp"
#include "Omega_h_string.hpp"

namespace Omega_h {

namespace {

void print_indicator(
    std::ostream& os, std::string const& above, std::size_t pos) {
  for (std::size_t i = 0; i < pos; ++i) {
    if (above.at(i) == '\t')
      os << '\t';
    else
      os << ' ';
  }
  os << "^\n";
}

void print_underline(std::ostream& os, std::string const& above,
    std::size_t start, std::size_t end) {
  for (std::size_t i = 0; i < start; ++i) {
    if (above.at(i) == '\t')
      os << '\t';
    else
      os << ' ';
  }
  for (auto i = start; i < end; ++i) os << '~';
  os << '\n';
}

}  // end anonymous namespace

void Reader::at_token(std::istream& stream) {
  bool done = false;
  /* this can loop arbitrarily as reductions are made,
     because they don't consume the token */
  while (!done) {
    auto parser_action = get_action(parser, parser_state, lexer_token);
    if (parser_action.kind == ACTION_NONE) {
      std::stringstream ss;
      ss << "error: Parser failure at line " << line;
      ss << " column " << column << " of " << stream_name << '\n';
      error_print_line(stream, ss);
      std::set<std::string> expect_names;
      for (int expect_token = 0; expect_token < grammar->nterminals;
           ++expect_token) {
        auto expect_action = get_action(parser, parser_state, expect_token);
        if (expect_action.kind != ACTION_NONE) {
          expect_names.insert(at(grammar->symbol_names, expect_token));
        }
      }
      ss << "Expected one of {";
      for (std::set<std::string>::iterator it = expect_names.begin();
           it != expect_names.end(); ++it) {
        if (it != expect_names.begin()) ss << ", ";
        if (*it == ",")
          ss << "','";
        else
          ss << *it;
      }
      ss << "}\n";
      ss << "Got: " << at(grammar->symbol_names, lexer_token) << '\n';
      ss << "Lexer text: \"" << lexer_text << "\"\n";
      ss << "Parser was in state " << parser_state << '\n';
      throw ParserFail(ss.str());
    } else if (parser_action.kind == ACTION_SHIFT) {
      if (sensing_indent) {
        symbol_indentation_stack.push_back(indent_text.size());
      }
      Omega_h::any shift_result;
      shift_result = this->at_shift(lexer_token, lexer_text);
      value_stack.emplace_back(std::move(shift_result));
      done = true;
    } else if (parser_action.kind == ACTION_REDUCE) {
      if (parser_action.production == get_accept_production(*grammar)) {
        did_accept = true;
        return;
      }
      auto& prod = at(grammar->productions, parser_action.production);
      reduction_rhs.clear();
      for (int i = 0; i < size(prod.rhs); ++i) {
        reduction_rhs.emplace_back(
            std::move(at(value_stack, size(value_stack) - size(prod.rhs) + i)));
      }
      resize(value_stack, size(value_stack) - size(prod.rhs));
      Omega_h::any reduce_result;
      try {
        reduce_result =
            this->at_reduce(parser_action.production, reduction_rhs);
      } catch (const ParserFail& e) {
        std::stringstream ss;
        ss << "error: Parser failure at line " << line;
        ss << " column " << column << " of " << stream_name << '\n';
        error_print_line(stream, ss);
        ss << '\n' << e.what();
        throw ParserFail(ss.str());
      }
      value_stack.emplace_back(std::move(reduce_result));
      if (sensing_indent) {
        if (size(prod.rhs)) {
          resize(symbol_indentation_stack,
              (size(symbol_indentation_stack) + 1) - size(prod.rhs));
        } else {
          symbol_indentation_stack.push_back(symbol_indentation_stack.back());
        }
      }
    } else {
      Omega_h_fail("SERIOUS BUG: Action::kind enum value not in range\n");
    }
    parser_state = execute_action(parser, parser_stack, parser_action);
  }
}

void Reader::indent_mismatch() {
  OMEGA_H_CHECK(!indent_stack.empty());
  auto top = indent_stack.back();
  std::stringstream ss;
  ss << "error: Indentation characters beginning line " << line << " of "
     << stream_name << " don't match those beginning line " << top.line << '\n';
  ss << "It is strongly recommended not to mix tabs and spaces in "
        "indentation-sensitive formats\n";
  throw ParserFail(ss.str());
}

void Reader::at_token_indent(std::istream& stream) {
  if (!sensing_indent || lexer_token != tables->indent_info.newline_token) {
    at_token(stream);
    return;
  }
  auto last_newline_pos = lexer_text.find_last_of("\n");
  if (last_newline_pos == std::string::npos) {
    Omega_h_fail("INDENT token did not contain a newline '\\n' !\n");
  }
  auto lexer_indent =
      lexer_text.substr(last_newline_pos + 1, std::string::npos);
  // the at_token call is allowed to do anything to lexer_text
  at_token(stream);
  lexer_text.clear();
  std::size_t minlen = std::min(lexer_indent.length(), indent_text.length());
  if (lexer_indent.length() > indent_text.length()) {
    if (0 != lexer_indent.compare(0, indent_text.length(), indent_text)) {
      indent_mismatch();
    }
    indent_stack.push_back({line, indent_text.length(), lexer_indent.length()});
    indent_text = lexer_indent;
    lexer_token = tables->indent_info.indent_token;
    at_token(stream);
  } else if (lexer_indent.length() < indent_text.length()) {
    if (0 != indent_text.compare(0, lexer_indent.length(), lexer_indent)) {
      indent_mismatch();
    }
    while (!indent_stack.empty()) {
      auto top = indent_stack.back();
      if (top.end_length <= minlen) break;
      indent_stack.pop_back();
      lexer_token = tables->indent_info.dedent_token;
      at_token(stream);
    }
    indent_text = lexer_indent;
  } else {
    if (0 != lexer_indent.compare(indent_text)) {
      indent_mismatch();
    }
  }
}

void Reader::backtrack_to_last_accept(std::istream& stream) {
  /* all the last_accept and backtracking is driven by
    the "accept the longest match" rule */
  line = last_lexer_accept_line;
  column = last_lexer_accept_column;
  line_text = last_lexer_accept_line_text;
  while (lexer_text.size() > last_lexer_accept) {
    bool ok = !stream.unget().fail();
    OMEGA_H_CHECK(ok);
    lexer_text.pop_back();
  }
}

void Reader::reset_lexer_state() {
  lexer_state = 0;
  lexer_text.clear();
  lexer_token = -1;
}

void Reader::at_lexer_end(std::istream& stream) {
  if (lexer_token == -1) {
    std::stringstream ss;
    if (lexer_text.find('\n') == std::string::npos) {
      ss << "error: Could not tokenize this (line " << line;
      ss << " column " << column << " of " << stream_name << "):\n";
      ss << line_text << '\n';
      OMEGA_H_CHECK(line_text.size() >= lexer_text.size());
      print_underline(ss, line_text, line_text.size() - lexer_text.size(),
          line_text.size());
    } else {
      ss << "error: Could not tokenize this (ends at line " << line;
      ss << " column " << column << " of " << stream_name << "):\n";
      ss << lexer_text << '\n';
    }
    throw ParserFail(ss.str());
  }
  backtrack_to_last_accept(stream);
  at_token_indent(stream);
  reset_lexer_state();
}

Reader::Reader(ReaderTablesPtr tables_in)
    : tables(tables_in),
      parser(tables->parser),
      lexer(tables->lexer),
      grammar(get_grammar(parser)) {
  OMEGA_H_CHECK(get_determinism(lexer));
}

void Reader::update_position(char c) {
  if (c == '\n') {
    ++line;
    column = 1;
    line_text.clear();
  } else {
    ++column;
  }
}

void Reader::error_print_line(std::istream& is, std::ostream& os) {
  auto oldpos = line_text.size();
  char c;
  while (is.get(c)) {
    if (c == '\n' || c == '\r') break;
    line_text.push_back(c);
  }
  if (line_text.empty()) return;
  os << line_text << '\n';
  if (oldpos > 0) print_indicator(os, line_text, oldpos - 1);
}

any Reader::read_stream(
    std::istream& stream, std::string const& stream_name_in) {
  line = 1;
  column = 1;
  lexer_state = 0;
  lexer_text.clear();
  line_text.clear();
  lexer_token = -1;
  parser_state = 0;
  parser_stack.clear();
  parser_stack.push_back(parser_state);
  value_stack.clear();
  did_accept = false;
  stream_name = stream_name_in;
  if (tables->indent_info.is_sensitive) {
    sensing_indent = true;
    indent_text.clear();
    indent_stack.clear();
  } else {
    sensing_indent = false;
  }
  char c;
  while (stream.get(c)) {
    if (!is_symbol(c)) {
      std::stringstream ss;
      ss << "error: Unexpected character code " << int(c);
      ss << " at line " << line << " column " << column;
      ss << " of " << stream_name << '\n';
      error_print_line(stream, ss);
      throw ParserFail(ss.str());
    }
    line_text.push_back(c);
    lexer_text.push_back(c);
    auto lexer_symbol = get_symbol(c);
    lexer_state = step(lexer, lexer_state, lexer_symbol);
    if (lexer_state == -1) {
      at_lexer_end(stream);
    } else {
      auto token = accepts(lexer, lexer_state);
      update_position(c);
      if (token != -1) {
        lexer_token = token;
        last_lexer_accept = lexer_text.size();
        last_lexer_accept_line = line;
        last_lexer_accept_column = column;
        last_lexer_accept_line_text = line_text;
      }
    }
  }
  if (last_lexer_accept < lexer_text.size()) {
    std::stringstream ss;
    std::string bad_str =
        lexer_text.substr(last_lexer_accept, std::string::npos);
    ss << "error: Could not tokenize \"" << bad_str;
    ss << "\" at end of " << stream_name << '\n';
    throw ParserFail(ss.str());
  }
  at_lexer_end(stream);
  lexer_token = get_end_terminal(*grammar);
  at_token(stream);
  if (!did_accept) {
    Omega_h_fail(
        "The EOF terminal was accepted but the root nonterminal was not "
        "reduced\n"
        "This indicates a bug in Omega_h::Reader\n");
  }
  OMEGA_H_CHECK(value_stack.size() == 1);
  return std::move(value_stack.back());
}

any Reader::read_string(
    std::string const& string, std::string const& string_name) {
  std::istringstream stream(string);
  return read_stream(stream, string_name);
}

any Reader::read_file(std::string const& file_name) {
  std::ifstream stream(file_name.c_str());
  if (!stream.is_open()) {
    std::stringstream ss;
    ss << "Could not open file " << file_name;
    throw ParserFail(ss.str());
  }
  return read_stream(stream, file_name);
}

any Reader::at_shift(int, std::string&) { return any(); }

any Reader::at_reduce(int, std::vector<any>&) { return any(); }

DebugReader::DebugReader(ReaderTablesPtr tables_in, std::ostream& os_in)
    : Reader(tables_in), os(os_in) {}

any DebugReader::at_shift(int token, std::string& text) {
  std::string text_escaped;
  for (auto c : text) {
    switch (c) {
      case '\n':
        text_escaped.append("\\n");
        break;
      case '\t':
        text_escaped.append("\\t");
        break;
      case '\r':
        text_escaped.append("\\r");
        break;
      default:
        text_escaped.push_back(c);
    }
  }
  os << "SHIFT (" << at(grammar->symbol_names, token) << ")[" << text_escaped
     << "]\n";
  return any(std::move(text_escaped));
}

any DebugReader::at_reduce(int prod_i, std::vector<any>& rhs) {
  os << "REDUCE";
  std::string lhs_text;
  auto& prod = at(grammar->productions, prod_i);
  for (int i = 0; i < size(prod.rhs); ++i) {
    auto& rhs_name = at(grammar->symbol_names, at(prod.rhs, i));
    auto rhs_text = move_value<std::string>(at(rhs, i));
    os << " (" << rhs_name << ")[" << rhs_text << "]";
    lhs_text.append(rhs_text);
  }
  auto& lhs_name = at(grammar->symbol_names, prod.lhs);
  os << " -> (" << lhs_name << ")[" << lhs_text << "]\n";
  return any(std::move(lhs_text));
}

}  // end namespace Omega_h
