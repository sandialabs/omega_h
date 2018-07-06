#ifndef OMEGA_H_EXPR_HPP
#define OMEGA_H_EXPR_HPP

#include <functional>
#include <map>
#include <vector>

#include <Omega_h_array.hpp>
#include <Omega_h_matrix.hpp>

#include <Omega_h_teuchos_includes.hpp>

namespace Omega_h {

class ExprReader : public Teuchos::Reader {
 public:
  using Args = std::vector<Teuchos::any>;
  using Function = std::function<void(Teuchos::any&, Args&)>;

 private:
  LO size;
  Int dim;
  std::map<std::string, Teuchos::any> variables;
  std::map<std::string, Function> functions;

 public:
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
};

}  // end namespace Omega_h

#endif
