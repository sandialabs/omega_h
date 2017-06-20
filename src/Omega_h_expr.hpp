#ifndef OMEA_H_EXPR_HPP
#define OMEA_H_EXPR_HPP

#include <functional>
#include <map>
#include <vector>

#include <Omega_h_array.hpp>
#include <Teuchos_Reader.hpp>

namespace Omega_h {


class ExprReader : public Teuchos::Reader {
 public:
  using Args = std::vector<Teuchos::any>;
  using Function = std::function<void(Teuchos::any&, Args&)>;
 private:
  LO count;
  Int dim;
  std::map<std::string, Teuchos::any> vars;
  std::map<std::string, Function> functions;
 public:
  ExprReader(LO count_in, Int dim_in);
  virtual ~ExprReader() override final;
  void register_variable(std::string const& name, Teuchos::any& value);
  void register_function(std::string const& name, Function const& value);
 protected:
  void at_shift(Teuchos::any& result, int token, std::string& text) override final;
  void at_reduce(Teuchos::any& result, int token, std::string& text) override final;
};

}  // end namespace Omega_h

#endif
