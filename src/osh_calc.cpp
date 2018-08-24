#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <Omega_h_language.hpp>
#include <Omega_h_reader.hpp>
#include <Omega_h_fail.hpp>

namespace {

enum {
  PROD_EXPR,
  PROD_ADD_SUB_DECAY,
  PROD_MUL_DIV_DECAY,
  PROD_POW_DECAY,
  PROD_NEG_DECAY,
  PROD_ADD,
  PROD_SUB,
  PROD_MUL,
  PROD_DIV,
  PROD_POW,
  PROD_UNARY_CALL,
  PROD_BINARY_CALL,
  PROD_PARENS,
  PROD_CONST,
  PROD_NEG,
  PROD_NO_SPACES,
  PROD_SPACES,
};

enum { NPRODS = PROD_SPACES + 1 };

enum {
  TOK_SPACE,
  TOK_NAME,
  TOK_ADD,
  TOK_SUB,
  TOK_MUL,
  TOK_DIV,
  TOK_POW,
  TOK_LPAREN,
  TOK_RPAREN,
  TOK_COMMA,
  TOK_CONST,
};

enum { NTOKS = TOK_CONST + 1 };

Omega_h::Language build_language() {
  Omega_h::Language out;
  auto& prods = out.productions;
  prods.resize(NPRODS);
  prods[PROD_EXPR] = {"expr", {"expr(+-)"}};
  prods[PROD_ADD_SUB_DECAY] = {"expr(+-)", {"expr(*/)"}};
  prods[PROD_MUL_DIV_DECAY] = {"expr(*/)", {"expr(^)"}};
  prods[PROD_POW_DECAY] = {"expr(^)", {"neg-expr"}};
  prods[PROD_NEG_DECAY] = {"neg-expr", {"scalar-expr"}};
  prods[PROD_ADD] = {"expr(+-)", {"expr(+-)", "+", "S?", "expr(*/)"}};
  prods[PROD_SUB] = {"expr(+-)", {"expr(+-)", "-", "S?", "expr(*/)"}};
  prods[PROD_MUL] = {"expr(*/)", {"expr(*/)", "*", "S?", "expr(^)"}};
  prods[PROD_DIV] = {"expr(*/)", {"expr(*/)", "/", "S?", "expr(^)"}};
  prods[PROD_POW] = {"expr(^)", {"expr(^)", "^", "S?", "neg-expr"}};
  prods[PROD_NEG] = {"neg-expr", {"-", "scalar-expr"}};
  prods[PROD_UNARY_CALL] = {"scalar-expr", {"name", "S?", "(", "S?", "expr(+-)", ")", "S?"}};
  prods[PROD_BINARY_CALL] = {"scalar-expr", {"name", "S?", "(", "S?", "expr(+-)", ",", "S?", "expr(+-)", ")", "S?"}};
  prods[PROD_PARENS] = {"scalar-expr", {"(", "S?", "expr(+-)", ")", "S?"}};
  prods[PROD_CONST] = {"scalar-expr", {"constant", "S?"}};
  prods[PROD_NO_SPACES] = {"S?", {}};
  prods[PROD_SPACES] = {"S?", {"spaces"}};
  out.tokens.resize(NTOKS);
  out.tokens[TOK_SPACE] = {"spaces", "[ \t\n\r]+"};
  out.tokens[TOK_NAME] = {"name", "[_a-zA-Z][_a-zA-Z0-9]*"};
  out.tokens[TOK_ADD] = {"+", "\\+"};
  out.tokens[TOK_SUB] = {"-", "\\-"};
  out.tokens[TOK_MUL] = {"*", "\\*"};
  out.tokens[TOK_DIV] = {"/", "\\/"};
  out.tokens[TOK_POW] = {"^", "\\^"};
  out.tokens[TOK_LPAREN] = {"(", "\\("};
  out.tokens[TOK_RPAREN] = {")", "\\)"};
  out.tokens[TOK_COMMA] = {",", ","};
  out.tokens[TOK_CONST] = {"constant", "(0|([1-9][0-9]*))(\\.[0-9]*)?([eE]\\-?[1-9][0-9]*)?"};
  return out;
}

Omega_h::ReaderTablesPtr ask_reader_tables() {
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif
  static Omega_h::ReaderTablesPtr ptr;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
  if (ptr.use_count() == 0) {
    ptr = Omega_h::build_reader_tables(build_language());
  }
  return ptr;
}

class Reader : public Omega_h::Reader {
 public:
  Reader():Omega_h::Reader(ask_reader_tables()) {
    unary_map["sqrt"] = &std::sqrt;
    unary_map["sin"] = &std::sin;
    unary_map["cos"] = &std::cos;
    unary_map["tan"] = &std::tan;
    unary_map["asin"] = &std::asin;
    unary_map["acos"] = &std::acos;
    unary_map["atan"] = &std::atan;
    unary_map["exp"] = &std::exp;
    unary_map["log"] = &std::log;
    unary_map["log10"] = &std::log10;
    unary_map["log2"] = &std::log2;
    unary_map["cbrt"] = &std::cbrt;
    binary_map["hypot"] = &std::hypot;
    binary_map["atan2"] = &std::atan2;
  }
  Reader(Reader const& other) = default;
  virtual ~Reader() override = default;
 protected:
  virtual Omega_h::any at_shift(int token, std::string& text) override {
    switch (token) {
      case TOK_NAME: return Omega_h::any(std::move(text));
      case TOK_CONST: {
        return Omega_h::any(std::atof(text.c_str()));
      }
    }
    return Omega_h::any();
  }
  virtual Omega_h::any at_reduce(int prod, std::vector<Omega_h::any>& rhs) override {
    switch (prod) {
      case PROD_EXPR:
      case PROD_ADD_SUB_DECAY:
      case PROD_MUL_DIV_DECAY:
      case PROD_POW_DECAY:
      case PROD_NEG_DECAY:
        return std::move(rhs.at(0));
      case PROD_ADD:
        return Omega_h::move_value<double>(rhs.at(0)) + Omega_h::move_value<double>(rhs.at(3));
      case PROD_SUB:
        return Omega_h::move_value<double>(rhs.at(0)) - Omega_h::move_value<double>(rhs.at(3));
      case PROD_MUL:
        return Omega_h::move_value<double>(rhs.at(0)) * Omega_h::move_value<double>(rhs.at(3));
      case PROD_DIV:
        return Omega_h::move_value<double>(rhs.at(0)) / Omega_h::move_value<double>(rhs.at(3));
      case PROD_POW:
        return std::pow(Omega_h::move_value<double>(rhs.at(0)), Omega_h::move_value<double>(rhs.at(3)));
      case PROD_NEG:
        return - Omega_h::move_value<double>(rhs.at(1));
      case PROD_UNARY_CALL: {
        auto name = Omega_h::move_value<std::string>(rhs.at(0));
        auto arg = Omega_h::move_value<double>(rhs.at(4));
        if (!unary_map.count(name)) {
          std::stringstream ss;
          ss << "Unknown unary function name \"" << name << "\"\n";
          throw Omega_h::ParserFail(ss.str());
        }
        auto fptr = unary_map[name];
        return (*fptr)(arg);
      }
      case PROD_BINARY_CALL: {
        auto name = Omega_h::move_value<std::string>(rhs.at(0));
        auto arg1 = Omega_h::move_value<double>(rhs.at(4));
        auto arg2 = Omega_h::move_value<double>(rhs.at(7));
        if (!binary_map.count(name)) {
          std::stringstream ss;
          ss << "Unknown binary function name \"" << name << "\"\n";
          throw Omega_h::ParserFail(ss.str());
        }
        auto fptr = binary_map[name];
        return (*fptr)(arg1, arg2);
      }
      case PROD_PARENS:
        return std::move(rhs.at(2));
      case PROD_CONST:
        return std::move(rhs.at(0));
    }
    OMEGA_H_NORETURN(Omega_h::any());
  }
 private:
  typedef double (*Unary)(double);
  typedef double (*Binary)(double, double);
  std::map<std::string, Unary> unary_map;
  std::map<std::string, Binary> binary_map;
};

}  // end anonymous namespace

int main() {
  auto reader = Reader();
  for (std::string line; std::getline(std::cin, line);) {
    reader.read_string("input", line);
    auto value = reader.move_result<double>();
    std::cout << value << '\n';
  }
}
