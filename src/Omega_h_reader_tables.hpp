#ifndef OMEGA_H_READER_TABLES_HPP
#define OMEGA_H_READER_TABLES_HPP

#include <memory>

#include <Omega_h_finite_automaton.hpp>
#include <Omega_h_parser.hpp>

namespace Omega_h {

struct IndentInfo {
  bool is_sensitive;
  int indent_token;
  int dedent_token;
  int newline_token;
};

struct ReaderTables {
  Parser parser;
  FiniteAutomaton lexer;
  IndentInfo indent_info;
};

using ReaderTablesPtr = std::shared_ptr<ReaderTables const>;

}  // namespace Omega_h

#endif
