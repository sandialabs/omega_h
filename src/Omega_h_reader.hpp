#ifndef OMEGA_H_READER_HPP
#define OMEGA_H_READER_HPP

#include <functional>
#include <iosfwd>

#include <Omega_h_any.hpp>
#include <Omega_h_reader_tables.hpp>
#include <Omega_h_std_vector.hpp>

namespace Omega_h {

class Reader {
 public:
  Reader() = delete;
  Reader(Reader const&) = default;
  virtual ~Reader() = default;
  Reader(ReaderTablesPtr tables_in);
  any read_stream(std::istream& stream, std::string const& stream_name_in = "");
  any read_string(
      std::string const& string, std::string const& string_name = "");
  any read_file(std::string const& file_name);

 protected:
  virtual any at_shift(int token, std::string& text);
  virtual any at_reduce(int token, std::vector<any>& rhs);

 protected:
  ReaderTablesPtr tables;
  Parser const& parser;
  FiniteAutomaton const& lexer;
  GrammarPtr grammar;
  std::size_t line;
  std::size_t column;
  int lexer_state;
  std::string lexer_text;
  std::string line_text;
  int lexer_token;
  std::size_t last_lexer_accept;
  std::size_t last_lexer_accept_line;
  std::size_t last_lexer_accept_column;
  std::string last_lexer_accept_line_text;
  int parser_state;
  std::vector<int> parser_stack;
  std::vector<any> value_stack;
  std::vector<any> reduction_rhs;
  std::string stream_name;
  bool did_accept;

 protected:  // variables for indentation-sensitive language parsing
  bool sensing_indent;
  std::string indent_text;
  struct IndentStackEntry {
    std::size_t line;
    std::size_t start_length;
    std::size_t end_length;
  };
  // this is the stack that shows, for the current leading indentation
  // characters, which subset of them came from each nested increase
  // in indentation
  std::vector<IndentStackEntry> indent_stack;
  // this stack notes, for each symbol in the pushdown automaton
  // stack, how many characters indent the line that that symbol
  // starts on
  std::vector<std::size_t> symbol_indentation_stack;

 private:  // helper methods
  void at_token(std::istream& stream);
  [[noreturn]] void indent_mismatch();
  void at_token_indent(std::istream& stream);
  void at_lexer_end(std::istream& stream);
  void backtrack_to_last_accept(std::istream& stream);
  void reset_lexer_state();
  void update_position(char c);
  void error_print_line(std::istream& is, std::ostream& os);
};

class DebugReader : public Reader {
 public:
  DebugReader(ReaderTablesPtr tables_in, std::ostream& os_in);
  DebugReader(DebugReader const& other) = default;
  virtual ~DebugReader() override = default;

 protected:
  virtual any at_shift(int token, std::string& text) override;
  virtual any at_reduce(int token, std::vector<any>& rhs) override;

 private:
  std::ostream& os;
};

}  // namespace Omega_h

#endif
