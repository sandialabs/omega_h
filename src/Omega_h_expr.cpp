#include <Teuchos_Language.hpp>
#include <Teuchos_MathExpr.hpp>

#include <Omega_h_array_ops.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_vector.hpp>

namespace Omega_h {

using Teuchos::any;
using Teuchos::any_cast;

namespace {

void promote_bool(LO size, any& x) {
  x = Bytes(size, Byte(any_cast<bool>(x)));
}

void promote_bools(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(bool) && rhs.type() == typeid(Bytes)) {
    TEUCHOS_ASSERT(any_cast<Bytes>(rhs).size() == size);
    promote_bool(size, lhs);
  }
  if (lhs.type() == typeid(Bytes) && rhs.type() == typeid(bool)) {
    TEUCHOS_ASSERT(any_cast<Bytes>(lhs).size() == size);
    promote_bool(size, rhs);
  }
}

void promote_scalar(LO size, any& x) {
  x = Reals(size, any_cast<Real>(x));
}

void promote_scalars(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real) && rhs.type() == typeid(Reals)) {
    promote_scalar(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Real)) {
    promote_scalar(size, rhs);
  }
}

template <Int dim>
void promote_vector(LO size, any& x) {
  x = repeat_vector(size, any_cast<Vector<dim>>(x));
}

template <Int dim>
void promote_vectors(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Reals)) {
    promote_vector<dim>(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Vector<dim>)) {
    promote_vector<dim>(size, rhs);
  }
}

template <Int dim>
void promote_matrix(LO size, any& x) {
  x = repeat_matrix(size, any_cast<Matrix<dim,dim>>(x));
}

template <Int dim>
void promote_matrices(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Matrix<dim, dim>) && rhs.type() == typeid(Reals)) {
    promote_matrix<dim>(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Matrix<dim, dim>)) {
    promote_matrix<dim>(size, rhs);
  }
}

template <Int dim>
void promote(LO size, any& lhs, any& rhs) {
  if (lhs.type() == rhs.type()) return;
  promote_bools(size, lhs, rhs);
  promote_scalars(size, lhs, rhs);
  promote_vectors<dim>(size, lhs, rhs);
  promote_matrices<dim>(size, lhs, rhs);
}

void promote(LO size, Int dim, any& lhs, any& rhs) {
  if (dim == 3) promote<3>(size, lhs, rhs);
  else if (dim == 2) promote<2>(size, lhs, rhs);
  else if (dim == 1) promote<1>(size, lhs, rhs);
  OMEGA_H_NORETURN();
}

template <Int dim>
void promote(LO size, any& x) {
  promote_bool(size, x);
  promote_scalar(size, x);
  promote_vector<dim>(size, x);
  promote_matrix<dim>(size, x);
}

void promote(LO size, Int dim, any& x) {
  if (dim == 3) promote<3>(size, x);
  else if (dim == 2) promote<2>(size, x);
  else if (dim == 1) promote<1>(size, x);
  OMEGA_H_NORETURN();
}

any ternary(LO size, Int dim, any& cond,
    any& true_val, any& false_val) {
  if (cond.type() == typeid(bool)) {
    if (true_val.type() == typeid(Real)) {
      return any_cast<bool>(cond) ?
             any_cast<Real>(true_val) :
             any_cast<Real>(false_val);
    } else if (true_val.type() == typeid(Reals)) {
      return any_cast<bool>(cond) ?
             any_cast<Reals>(true_val) :
             any_cast<Reals>(false_val);
    } else {
      throw Teuchos::ParserFail("Invalid true value type in ternary operator");
    }
  } else if (cond.type() == typeid(Bytes)) {
    promote(size, dim, true_val);
    promote(size, dim, false_val);
    return Reals(ternary_each(any_cast<Bytes>(cond),
               any_cast<Reals>(true_val),
               any_cast<Reals>(false_val)));
  } else {
    throw Teuchos::ParserFail("Invalid condition value type in ternary operator");
  }
}

any eval_or(any& lhs, any& rhs) {
  if (lhs.type() == typeid(bool)) {
    return any_cast<bool>(lhs) || any_cast<bool>(rhs);
  } else if (lhs.type() == typeid(Bytes)) {
    return lor_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to || operator");
  }
}

any eval_and(any& lhs, any& rhs) {
  if (lhs.type() == typeid(bool)) {
    return any_cast<bool>(lhs) && any_cast<bool>(rhs);
  } else if (lhs.type() == typeid(Bytes)) {
    return land_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to && operator");
  }
}

any gt(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) > any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = gt_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to > operator");
  }
}

any lt(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) < any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = lt_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to > operator");
  }
}

any eq(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) == any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = eq_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to == operator");
  }
}

template <Int dim>
any add(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) + any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Vector<dim>>(lhs) + any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Matrix<dim,dim>)) {
    result = any_cast<Matrix<dim,dim>>(lhs) + any_cast<Matrix<dim,dim>>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = Reals(add_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to + operator");
  }
}

any add(Int dim, any& lhs, any& rhs) {
  if (dim == 3) return add<3>(lhs, rhs);
  if (dim == 2) return add<2>(lhs, rhs);
  if (dim == 1) return add<1>(lhs, rhs);
  OMEGA_H_NORETURN(any);
}

template <Int dim>
any sub(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) - any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Vector<dim>>(lhs) - any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Matrix<dim,dim>)) {
    result = any_cast<Matrix<dim,dim>>(lhs) - any_cast<Matrix<dim,dim>>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = Reals(subtract_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to + operator");
  }
}

any sub(Int dim, any& lhs, any& rhs) {
  if (dim == 3) return sub<3>(lhs, rhs);
  if (dim == 2) return sub<2>(lhs, rhs);
  if (dim == 1) return sub<1>(lhs, rhs);
  OMEGA_H_NORETURN(any);
}

template <Int dim>
any mul(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real) && rhs.type() == typeid(Real)) {
    return any_cast<Real>(lhs) * any_cast<Real>(rhs);
  /* begin multiply non-scalar by scalar (commutative) */
  } else if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Real)) {
    return any_cast<Vector<dim>>(lhs) * any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Real) && rhs.type() == typeid(Vector<dim>)) {
    return any_cast<Real>(lhs) * any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Matrix<dim,dim>) && rhs.type() == typeid(Real)) {
    return any_cast<Matrix<dim,dim>>(lhs) * any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Real) && rhs.type() == typeid(Matrix<dim,dim>)) {
    return any_cast<Real>(lhs) * any_cast<Matrix<dim,dim>>(rhs);
  /* dot product */
  } else if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Vector<dim>)) {
    return any_cast<Vector<dim>>(lhs) * any_cast<Vector<dim>>(rhs);
  /* matrix * vector (non-commutative) */
  } else if (lhs.type() == typeid(Matrix<dim,dim>) && rhs.type() == typeid(Vector<dim>)) {
    return any_cast<Matrix<dim,dim>>(lhs) * any_cast<Vector<dim>>(rhs);
  /* matrix * matrix (non-commutative) */
  } else if (lhs.type() == typeid(Matrix<dim,dim>) &&
      rhs.type() == typeid(Matrix<dim,dim>)) {
    return any_cast<Matrix<dim,dim>>(lhs) * any_cast<Matrix<dim,dim>>(rhs);
  } else if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Reals)) {
    auto& lhs_vals = any_cast<Reals>(lhs);
    auto& rhs_vals = any_cast<Reals>(rhs);
    if (rhs_vals.size() == size) { // RHS is scalars
      return Reals(multiply_each(lhs_vals, rhs_vals));
    } else if (lhs_vals.size() == size) { // LHS is scalars
      return Reals(multiply_each(lhs_vals, rhs_vals));
    } else if (lhs_vals.size() == size * dim &&
        rhs_vals.size() == size * dim) { // dot products
      return dot_vectors(lhs_vals, rhs_vals);
    } else if (lhs_vals.size() == size * matrix_ncomps(dim) &&
        rhs_vals.size() == size * dim) { // matrices * vectors
      return matrices_times_vectors(lhs_vals, rhs_vals, dim);
    } else if (lhs_vals.size() == size * matrix_ncomps(dim) &&
        rhs_vals.size() == size * matrix_ncomps(dim)) {
      return matrices_times_matrices(lhs_vals, rhs_vals, dim);
    } else {
      throw Teuchos::ParserFail("Unexpected array size in * operator");
    }
  } else {
    throw Teuchos::ParserFail("Invalid operand types to * operator");
  }
}

any mul(LO size, Int dim, any& lhs, any& rhs) {
  if (dim == 3) return mul<3>(size, lhs, rhs);
  if (dim == 2) return mul<2>(size, lhs, rhs);
  if (dim == 1) return mul<1>(size, lhs, rhs);
  OMEGA_H_NORETURN(any());
}

template <Int dim>
any div(LO size, any& lhs, any& rhs) {
  if (rhs.type() == typeid(Reals)) {
    return Reals(divide_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else if (rhs.type() == typeid(Real)) {
    if (lhs.type() == typeid(Real)) {
      return any_cast<Real>(lhs) / any_cast<Real>(rhs);
    } else if (lhs.type() == typeid(Vector<dim>)) {
      return any_cast<Vector<dim>>(lhs) / any_cast<Real>(rhs);
    } else if (lhs.type() == typeid(Matrix<dim,dim>)) {
      return any_cast<Matrix<dim,dim>>(lhs) / any_cast<Real>(rhs);
    } else {
      throw Teuchos::ParserFail("Invalid left operand type in / operator");
    }
  } else {
    throw Teuchos::ParserFail("Invalid right operand type in / operator");
  }
}

any div(LO size, Int dim, any& lhs, any& rhs) {
  if (dim == 3) return div<3>(size, lhs, rhs);
  if (dim == 2) return div<2>(size, lhs, rhs);
  if (dim == 1) return div<1>(size, lhs, rhs);
  OMEGA_H_NORETURN(any());
}

template <Int dim>
any eval_pow(LO size, any& lhs, any& rhs) {
  if (rhs.type() == typeid(Reals)) {
    return Reals(pow_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else if (rhs.type() == typeid(Real)) {
    return std::pow(any_cast<Real>(lhs), any_cast<Real>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid right operand type in ^ operator");
  }
}

any eval_pow(LO size, Int dim, any& lhs, any& rhs) {
  if (dim == 3) return eval_pow<3>(size, lhs, rhs);
  if (dim == 2) return eval_pow<2>(size, lhs, rhs);
  if (dim == 1) return eval_pow<1>(size, lhs, rhs);
  OMEGA_H_NORETURN(any());
}

template <Int dim>
any neg(any& val) {
  if (val.type() == typeid(Real)) {
    return - any_cast<Real>(val);
  } else if (val.type() == typeid(Vector<dim>)) {
    return - any_cast<Vector<dim>>(val);
  } else if (val.type() == typeid(Matrix<dim,dim>)) {
    return - any_cast<Matrix<dim,dim>>(val);
  } else if (val.type() == typeid(Reals)) {
    return Reals(multiply_each_by(-1.0, any_cast<Reals>(val)));
  } else {
    throw Teuchos::ParserFail("Invalid operand type to negation operator");
  }
}

any neg(Int dim, any& val) {
  if (dim == 3) return neg<3>(val);
  if (dim == 2) return neg<2>(val);
  if (dim == 1) return neg<1>(val);
  OMEGA_H_NORETURN(any());
}

template <Int dim>
any access(LO size, any& var, Args& args) {
  auto i = static_cast<Int>(any_cast<Real>(args.at(0)));
  auto j = args.size() > 1 ? static_cast<Int>(any_cast<Real>(args.at(0))) : Int(-1);
  if (var.type() == typeid(Vector<dim>)) {
    return (any_cast<Vector<dim>>(var))(i);
  } else if (var.type() == typeid(Matrix<dim, dim>)) {
    return (any_cast<Matrix<dim,dim>>(var))(i,j);
  } else if (var.type() == typeid(Reals)) {
    auto array = any_cast<Reals>(var);
    if (array.size() == size * dim) {
      return Reals(get_component(array, dim, i));
    } else if (array.size() == size * matrix_ncomps(dim)) {
      return Reals(get_component(array, matrix_ncomps(dim), j * dim + i));
    } else {
      throw Teuchos::ParserFail("Unexpected array size in access operator\n");
    }
  } else {
    throw Teuchos::ParserFail("Unexpected variable type in access operator\n");
  }
}

any access(LO size, Int dim, any& var, Args& args) {
  if (dim == 3) return access<3>(size, var, args);
  if (dim == 2) return access<2>(size, var, args);
  if (dim == 1) return access<1>(size, var, args);
  OMEGA_H_NORETURN(any());
}

}  // end anonymous namespace

ExprReader::ExprReader(LO count_in, Int dim_in):
  Teuchos::Reader(Teuchos::MathExpr::ask_reader_tables()),
  count(count_in),
  dim(dim_in)
{
}

ExprReader::~ExprReader() {
}

void ExprReader::at_shift(any& result_any, int token, std::string& text) override final {
  using std::swap;
  switch (token) {
    case Teuchos::MathExpr::TOK_NAME {
      auto& result = Teuchos::make_any_ref<std::string>(result_any); 
      swap(result, text);
      break;
    }
    case Teuchos::MathExpr::TOK_CONST {
      result_any = std::atof(text.c_str());
    }
  }
}

void ExprReader::at_reduce(any& result, int token, std::string& text) override final {
  using std::swap;
  switch (prod) {
    case Teuchos::MathExpr::PROD_EXPR:
    case Teuchos::MathExpr::PROD_TERNARY_DECAY:
    case Teuchos::MathExpr::PROD_OR_DECAY:
    case Teuchos::MathExpr::PROD_AND_DECAY:
    case Teuchos::MathExpr::PROD_ADD_SUB_DECAY:
    case Teuchos::MathExpr::PROD_MUL_DIV_DECAY:
    case Teuchos::MathExpr::PROD_POW_DECAY:
    case Teuchos::MathExpr::PROD_NEG_DECAY:
    case Teuchos::MathExpr::PROD_SOME_ARGS:
    case Teuchos::MathExpr::PROD_CONST:
      swap(result, rhs.at(0));
      break;
    case Teuchos::MathExpr::PROD_TERNARY:
      promote(size, dim, rhs.at(3), rhs.at(6));
      result = ternary(size, dim, rhs.at(0), rhs.at(3), rhs.at(6));
      break;
    case Teuchos::MathExpr::PROD_OR:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = eval_or(rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_AND:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = eval_and(rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_GT:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = gt(rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_LT:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = lt(rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_GEQ:
    case Teuchos::MathExpr::PROD_LEQ:
      throw Teuchos::ParserFail("Operators <= and >= not supported yet");
      break;
    case Teuchos::MathExpr::PROD_EQ:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = eq(rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_ADD:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = add(dim, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_SUB:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = sub(dim, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_MUL:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = mul(size, dim, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_DIV:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = div(size, dim, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_POW:
      promote(size, dim, rhs.at(0), rhs.at(3));
      result = eval_pow(size, dim, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_CALL: {
      auto& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      auto& args = Teuchos::any_ref_cast<Args>(rhs.at(4));
      auto vit = vars.find(name);
      if (vit == vars.end()) {
        /* function call */
        auto fit = functions.find(name);
        TEUCHOS_TEST_FOR_EXCEPTION(fit == functions.end(), Teuchos::ParserFail,
            "\"" << name << "\" is neither a variable nor a function name\n");
        fit->second(result, args);
      } else {
        /* access operator for vector/matrix */
        auto& val = vit->second;
        result = access(size, dim, val, args);
      }
      break;
    }
    case Teuchos::MathExpr::PROD_NO_ARGS: {
      result = Args{};
      break;
    }
    case Teuchos::MathExpr::PROD_FIRST_ARG: {
      auto& args = Teuchos::make_any_ref<Args>(result);
      args.push_back(any());
      swap(args.back(), rhs.at(0));
      break;
    }
    case Teuchos::MathExpr::PROD_NEXT_ARG: {
      auto& args = Teuchos::any_ref_cast<Args>(rhs.at(0));
      args.push_back(any());
      swap(args.back(), rhs.at(3));
      swap(result, rhs.at(0));
      break;
    }
    case Teuchos::MathExpr::PROD_NEG:
      result = neg(dim, rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_PARENS:
      result = rhs.at(2);
      break;
    case Teuchos::MathExpr::PROD_VAR:
      auto& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      auto it = vars.find(name);
      TEUCHOS_TEST_FOR_EXCEPTION(it == vars.end(), Teuchos::ParserFail,
          "unknown variable name \"" << name << "\"\n");
      result = it->second;
      break;
  }
}

}  // end namespace Omega_h
