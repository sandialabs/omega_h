#ifndef OMEGA_H_BUILD_PARSER_HPP
#define OMEGA_H_BUILD_PARSER_HPP

#include <memory>
#include <set>

#include <Omega_h_parser.hpp>
#include <Omega_h_parser_graph.hpp>

namespace Omega_h {

struct Config {
  int production;
  int dot;
};

using Configs = std::vector<Config>;

using Context = std::set<int>;

/* nonterminal transitions will be stored as SHIFT
   actions while in progress */
struct ActionInProgress {
  Action action;
  Context context;
};

struct StateInProgress {
  std::vector<int> configs;
  std::vector<ActionInProgress> actions;
};

using StatesInProgress = std::vector<std::unique_ptr<StateInProgress>>;

struct StateConfig {
  int state;
  int config_in_state;
};

using StateConfigs = std::vector<StateConfig>;

struct ParserInProgress {
  StatesInProgress states;
  Configs configs;
  StateConfigs state_configs;
  ParserGraph states2state_configs;
  GrammarPtr grammar;
};

StateConfigs form_state_configs(StatesInProgress const& states);
ParserGraph form_states_to_state_configs(
    StateConfigs const& scs, StatesInProgress const& states);

void print_dot(std::string const& filepath, ParserInProgress const& pip);

ParserInProgress build_lalr1_parser(GrammarPtr grammar, bool verbose = false);

Parser accept_parser(ParserInProgress const& pip);

}  // namespace Omega_h

#endif
