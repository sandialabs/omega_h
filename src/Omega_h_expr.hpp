#ifndef OMEGA_H_EXPR_HPP
#define OMEGA_H_EXPR_HPP

#include <functional>
#include <map>
#include <vector>

#include <Omega_h_array.hpp>
#include <Omega_h_matrix.hpp>

#include <Omega_h_teuchos_includes.hpp>

namespace Omega_h {

struct ExprEnv {
  ExprEnv(LO size_in, Int dim_in);
  using Args = std::vector<Teuchos::any>;
  using Function = std::function<void(Teuchos::any&, Args&)>;
  void register_variable(std::string const& name, Teuchos::any const& value);
  void register_function(std::string const& name, Function const& value);
  void repeat(Teuchos::any& x);
  std::map<std::string, Teuchos::any> variables;
  std::map<std::string, Function> functions;
  LO size;
  Int dim;
};

struct ExprOp {
  virtual ~ExprOp();
  virtual void eval(ExprEnv& env, Teuchos::any& result) = 0;
};

class ExprOpsReader : public Teuchos::Reader {
 public:
  ExprOpsReader();
  virtual ~ExprOpsReader() override final = default;
  std::shared_ptr<ExprOp> read_ops(std::string const& str);

 protected:
  void at_shift(
      Teuchos::any& result, int token, std::string& text) override final;
  void at_reduce(Teuchos::any& result, int token,
      std::vector<Teuchos::any>& rhs) override final;
};

class ExprReader : public Teuchos::Reader {
 public:
  using Args = ExprEnv::Args;
  using Function = ExprEnv::Function;
  ExprReader(LO size_in, Int dim_in);
  virtual ~ExprReader() override final;
  void register_variable(std::string const& name, Teuchos::any const& value);
  void register_function(std::string const& name, Function const& value);
  void repeat(Teuchos::any& x);

 protected:
  void at_shift(
      Teuchos::any& result, int token, std::string& text) override final;
  void at_reduce(Teuchos::any& result, int token,
      std::vector<Teuchos::any>& rhs) override final;
  ExprEnv env;
};

}  // end namespace Omega_h

#endif
