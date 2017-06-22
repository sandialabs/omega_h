#ifndef OMEGA_H_EXPR_HPP
#define OMEGA_H_EXPR_HPP

#include <functional>
#include <map>
#include <vector>

#include <Teuchos_Reader.hpp>

#include <Omega_h_array.hpp>
#include <Omega_h_matrix.hpp>

#ifndef OMEGA_H_USE_TEUCHOSPARSER
#error "Can't use Omega_h_expr.hpp without Omega_h_USE_TeuchosParser=ON"
#endif

/* appease the non-standard crap in Teuchos::any */
namespace Teuchos {

template <Omega_h::Int dim>
inline
bool operator==(Omega_h::Vector<dim> const&, Omega_h::Vector<dim> const&) {
  return false;
}

template <Omega_h::Int dim>
inline
bool operator==(Omega_h::Matrix<dim,dim> const&, Omega_h::Matrix<dim,dim> const&) {
  return false;
}

inline
bool operator==(Omega_h::Reals const&, Omega_h::Reals const&) {
  return false;
}

inline
bool operator==(Omega_h::Bytes const&, Omega_h::Bytes const&) {
  return false;
}

inline
bool operator==(std::vector<Teuchos::any> const&, std::vector<Teuchos::any> const&) {
  return false;
}

template <Omega_h::Int dim>
inline
std::ostream& operator<<(std::ostream& os, Omega_h::Vector<dim> const&) {
  return os;
}

template <Omega_h::Int dim>
inline
std::ostream& operator<<(std::ostream& os, Omega_h::Matrix<dim,dim> const&) {
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, Omega_h::Reals const&) {
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, Omega_h::Bytes const&) {
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, std::vector<Teuchos::any> const&) {
  return os;
}

} /* done appeasing the non-standard crap in Teuchos::any */

namespace Omega_h {

class ExprReader : public Teuchos::Reader {
 public:
  using Args = std::vector<Teuchos::any>;
  using Function = std::function<void(Teuchos::any&, Args&)>;
 private:
  LO size;
  Int dim;
  std::map<std::string, Teuchos::any> vars;
  std::map<std::string, Function> functions;
 public:
  ExprReader(LO size_in, Int dim_in);
  virtual ~ExprReader() override final;
  void register_variable(std::string const& name, Teuchos::any& value);
  void register_function(std::string const& name, Function const& value);
 protected:
  void at_shift(Teuchos::any& result, int token, std::string& text) override final;
  void at_reduce(Teuchos::any& result, int token, std::vector<Teuchos::any>& rhs) override final;
};

}  // end namespace Omega_h

#endif
