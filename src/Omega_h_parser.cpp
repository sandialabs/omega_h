#include "Omega_h_parser.hpp"

namespace Omega_h {

Parser::Parser(GrammarPtr g, int nstates_reserve)
    : grammar(g),
      terminal_table(g->nterminals, nstates_reserve),
      nonterminal_table(get_nnonterminals(*g), nstates_reserve) {}

int get_nstates(Parser const& p) { return get_nrows(p.terminal_table); }

int add_state(Parser& p) {
  auto state = get_nstates(p);
  resize(p.terminal_table, state + 1, get_ncols(p.terminal_table));
  resize(p.nonterminal_table, state + 1, get_ncols(p.nonterminal_table));
  for (int t = 0; t < p.grammar->nterminals; ++t) {
    Action action;
    action.kind = ACTION_NONE;
    at(p.terminal_table, state, t) = action;
  }
  for (int nt = 0; nt < get_nnonterminals(*(p.grammar)); ++nt) {
    at(p.nonterminal_table, state, nt) = -1;
  }
  return state;
}

void add_terminal_action(Parser& p, int state, int terminal, Action action) {
  OMEGA_H_CHECK(at(p.terminal_table, state, terminal).kind == ACTION_NONE);
  OMEGA_H_CHECK(action.kind != ACTION_NONE);
  if (action.kind == ACTION_SHIFT) {
    OMEGA_H_CHECK(0 <= action.next_state);
    OMEGA_H_CHECK(action.next_state < get_nstates(p));
  } else {
    OMEGA_H_CHECK(0 <= action.production);
    OMEGA_H_CHECK(action.production < size(p.grammar->productions));
  }
  at(p.terminal_table, state, terminal) = action;
}

void add_nonterminal_action(
    Parser& p, int state, int nonterminal, int next_state) {
  OMEGA_H_CHECK(0 <= next_state);
  OMEGA_H_CHECK(next_state < get_nstates(p));
  OMEGA_H_CHECK(at(p.nonterminal_table, state, nonterminal) == -1);
  at(p.nonterminal_table, state, nonterminal) = next_state;
}

Action const& get_action(Parser const& p, int state, int terminal) {
  return at(p.terminal_table, state, terminal);
}

int execute_action(
    Parser const& p, std::vector<int>& stack, Action const& action) {
  OMEGA_H_CHECK(action.kind != ACTION_NONE);
  if (action.kind == ACTION_SHIFT) {
    stack.push_back(action.next_state);
  } else {
    auto& prod = at(p.grammar->productions, action.production);
    for (int i = 0; i < size(prod.rhs); ++i) stack.pop_back();
    OMEGA_H_CHECK(p.grammar.get());
    auto& grammar = *(p.grammar);
    auto nt = as_nonterminal(grammar, prod.lhs);
    OMEGA_H_CHECK(!stack.empty());
    auto next_state = at(p.nonterminal_table, stack.back(), nt);
    stack.push_back(next_state);
  }
  return stack.back();
}

GrammarPtr const& get_grammar(Parser const& p) { return p.grammar; }

ParserFail::ParserFail(const std::string& msg) : std::invalid_argument(msg) {}

void ParserFail::out_of_line_virtual_method() {}

}  // end namespace Omega_h
