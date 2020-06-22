#ifndef OMEGA_H_PARSER_HPP
#define OMEGA_H_PARSER_HPP

#include <stack>
#include <stdexcept>

#include <Omega_h_grammar.hpp>
#include <Omega_h_table.hpp>

namespace Omega_h {

enum ActionKind {
  ACTION_NONE,
  ACTION_SHIFT,
  ACTION_REDUCE,
};

struct Action {
  ActionKind kind;
  union {
    int production;
    int next_state;
  };
};

struct Parser {
  GrammarPtr grammar;
  /* (state x terminal) -> action */
  Table<Action> terminal_table;
  /* (state x non-terminal) -> new state */
  Table<int> nonterminal_table;
  Parser() = default;
  Parser(GrammarPtr g, int nstates_reserve);
};

int add_state(Parser& p);
int get_nstates(Parser const& p);
void add_terminal_action(Parser& p, int state, int terminal, Action action);
void add_nonterminal_action(
    Parser& p, int state, int nonterminal, int next_state);
Action const& get_action(Parser const& p, int state, int terminal);
int execute_action(
    Parser const& p, std::vector<int>& stack, Action const& action);
GrammarPtr const& get_grammar(Parser const& p);

class ParserFail : public std::invalid_argument {
 public:
  ParserFail(const std::string& msg);
  virtual void out_of_line_virtual_method();
};

}  // namespace Omega_h

#endif
