#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <Omega_h_fail.hpp>
#include <Omega_h_language.hpp>
#include <Omega_h_math_lang.hpp>
#include <Omega_h_reader.hpp>

namespace {

class CalcReader : public Omega_h::Reader {
 public:
  CalcReader() : Omega_h::Reader(Omega_h::math_lang::ask_reader_tables()) {
    unary_function_map["sqrt"] = &std::sqrt;
    unary_function_map["sin"] = &std::sin;
    unary_function_map["cos"] = &std::cos;
    unary_function_map["tan"] = &std::tan;
    unary_function_map["asin"] = &std::asin;
    unary_function_map["acos"] = &std::acos;
    unary_function_map["atan"] = &std::atan;
    unary_function_map["exp"] = &std::exp;
    unary_function_map["log"] = &std::log;
    unary_function_map["log10"] = &std::log10;
    binary_function_map["atan2"] = &std::atan2;
  }
  virtual ~CalcReader() = default;

 protected:
  struct CallArgs {
    double a0;
    double a1;
    int n;
  };
  virtual Omega_h::any at_shift(int token, std::string& text) {
    switch (token) {
      case Omega_h::math_lang::TOK_NAME: {
        return text;
      }
      case Omega_h::math_lang::TOK_CONST: {
        return std::atof(text.c_str());
      }
    }
    return Omega_h::any();
  }
  virtual Omega_h::any at_reduce(int prod, std::vector<Omega_h::any>& rhs) {
    using Omega_h::any_cast;
    switch (prod) {
      case Omega_h::math_lang::PROD_PROGRAM: {
        if (rhs.at(1).empty()) {
          throw Omega_h::ParserFail(
              "Calculator needs an expression to evaluate!");
        }
        return std::move(rhs.at(1));
      }
      case Omega_h::math_lang::PROD_NO_STATEMENTS:
      case Omega_h::math_lang::PROD_NO_EXPR:
      case Omega_h::math_lang::PROD_NEXT_STATEMENT: {
        return Omega_h::any();
      }
      case Omega_h::math_lang::PROD_ASSIGN: {
        auto& name = any_cast<std::string&>(rhs.at(0));
        double value = any_cast<double>(rhs.at(4));
        variable_map[name] = value;
        return value;
      }
      case Omega_h::math_lang::PROD_YES_EXPR:
      case Omega_h::math_lang::PROD_EXPR:
      case Omega_h::math_lang::PROD_TERNARY_DECAY:
      case Omega_h::math_lang::PROD_OR_DECAY:
      case Omega_h::math_lang::PROD_AND_DECAY:
      case Omega_h::math_lang::PROD_ADD_SUB_DECAY:
      case Omega_h::math_lang::PROD_MUL_DIV_DECAY:
      case Omega_h::math_lang::PROD_POW_DECAY:
      case Omega_h::math_lang::PROD_NEG_DECAY:
      case Omega_h::math_lang::PROD_SOME_ARGS:
        return std::move(rhs.at(0));
      case Omega_h::math_lang::PROD_TERNARY:
        return any_cast<bool>(rhs.at(0)) ? any_cast<double>(rhs.at(3))
                                         : any_cast<double>(rhs.at(6));
      case Omega_h::math_lang::PROD_OR:
        return any_cast<bool>(rhs.at(0)) || any_cast<bool>(rhs.at(3));
      case Omega_h::math_lang::PROD_AND:
        return any_cast<bool>(rhs.at(0)) && any_cast<bool>(rhs.at(3));
      case Omega_h::math_lang::PROD_GT:
        return any_cast<double>(rhs.at(0)) > any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_LT:
        return any_cast<double>(rhs.at(0)) < any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_GEQ:
        return any_cast<double>(rhs.at(0)) >= any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_LEQ:
        return any_cast<double>(rhs.at(0)) <= any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_EQ:
        return any_cast<double>(rhs.at(0)) == any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_BOOL_PARENS:
        return any_cast<bool>(rhs.at(2));
      case Omega_h::math_lang::PROD_ADD:
        return any_cast<double>(rhs.at(0)) + any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_SUB:
        return any_cast<double>(rhs.at(0)) - any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_MUL:
        return any_cast<double>(rhs.at(0)) * any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_DIV:
        return any_cast<double>(rhs.at(0)) / any_cast<double>(rhs.at(3));
      case Omega_h::math_lang::PROD_POW:
        return std::pow(
            any_cast<double>(rhs.at(0)), any_cast<double>(rhs.at(3)));
      case Omega_h::math_lang::PROD_CALL: {
        auto& name = any_cast<std::string&>(rhs.at(0));
        auto& args = any_cast<CallArgs&>(rhs.at(4));
        if (args.n < 1 || args.n > 2) {
          throw Omega_h::ParserFail(
              "Only unary and binary functions supported!\n");
        }
        if (args.n == 1) {
          if (!unary_function_map.count(name)) {
            std::stringstream ss;
            ss << "Unknown unary function name \"" << name << "\"\n";
            throw Omega_h::ParserFail(ss.str());
          }
          Unary fptr = unary_function_map[name];
          return (*fptr)(args.a0);
        } else {
          if (!binary_function_map.count(name)) {
            std::stringstream ss;
            ss << "Unknown binary function name \"" << name << "\"\n";
            throw Omega_h::ParserFail(ss.str());
          }
          Binary fptr = binary_function_map[name];
          return (*fptr)(args.a0, args.a1);
        }
      }
      case Omega_h::math_lang::PROD_NO_ARGS: {
        CallArgs args;
        args.n = 0;
        return args;
      }
      case Omega_h::math_lang::PROD_FIRST_ARG: {
        CallArgs args;
        args.a0 = any_cast<double>(rhs.at(0));
        args.n = 1;
        return args;
      }
      case Omega_h::math_lang::PROD_NEXT_ARG: {
        auto& args = any_cast<CallArgs&>(rhs.at(0));
        args.a1 = any_cast<double>(rhs.at(3));
        args.n = 2;
        return args;
      }
      case Omega_h::math_lang::PROD_NEG:
        return -any_cast<double>(rhs.at(2));
      case Omega_h::math_lang::PROD_VAL_PARENS:
        return any_cast<double>(rhs.at(2));
      case Omega_h::math_lang::PROD_CONST:
        return any_cast<double>(rhs.at(0));
      case Omega_h::math_lang::PROD_VAR:
        auto& name = any_cast<std::string&>(rhs.at(0));
        auto it = variable_map.find(name);
        if (it == variable_map.end()) {
          std::stringstream ss;
          ss << "variable " << name << " not defined!";
          throw Omega_h::ParserFail(ss.str());
        }
        return it->second;
    }
    return Omega_h::any();
  }

 private:
  typedef double (*Unary)(double);
  typedef double (*Binary)(double, double);
  std::map<std::string, Unary> unary_function_map;
  std::map<std::string, Binary> binary_function_map;
  std::map<std::string, double> variable_map;
};

}  // end anonymous namespace

int main() {
  CalcReader reader;
  for (std::string line; std::getline(std::cin, line);) {
    std::cout << Omega_h::any_cast<double>(reader.read_string("input", line))
              << '\n';
  }
}
