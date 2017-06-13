#include <Teuchos_Language.hpp>
#include <Teuchos_MathExpr.hpp>

namespace Omega_h {

using Teuchos::any;
using Teuchos::any_cast;

namespace {

void promote_bool(LO size, Teuchos::any& x) {
  x = Bytes(size, Byte(any_cast<bool>(x)));
}

void promote_bools(LO size, Teuchos::any& lhs, Teuchos::any& rhs) {
  if (lhs.type() == typeid(bool) && rhs.type() == typeid(Bytes)) {
    TEUCHOS_ASSERT(any_cast<Bytes>(rhs).size() == size);
    promote_bool(size, lhs);
  }
  if (lhs.type() == typeid(Bytes) && rhs.type() == typeid(bool)) {
    TEUCHOS_ASSERT(any_cast<Bytes>(lhs).size() == size);
    promote_bool(size, rhs);
  }
}

void promote_scalar(LO size, Teuchos::any& x) {
  x = Reals(size, any_cast<Real>(x));
}

void promote_scalars(LO size, Teuchos::any& lhs, Teuchos::any& rhs) {
  if (lhs.type() == typeid(Real) && rhs.type() == typeid(Reals)) {
    promote_scalar(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Real)) {
    promote_scalar(size, rhs);
  }
}

template <Int dim>
void promote_vector(LO size, Teuchos::any& x) {
  x = repeat_vector(size, any_cast<Vector<dim>>(x));
}

template <Int dim>
void promote_vectors(LO size, Teuchos::any& lhs, Teuchos::any& rhs) {
  if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Reals)) {
    promote_vector<dim>(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Vector<dim>)) {
    promote_vector<dim>(size, rhs);
  }
}

template <Int dim>
void promote_matrix(LO size, Teuchos::any& x) {
  x = repeat_matrix(size, any_cast<Matrix<dim,dim>>(x));
}

template <Int dim>
void promote_matrices(LO size, Teuchos::any& lhs, Teuchos::any& rhs) {
  if (lhs.type() == typeid(Matrix<dim, dim>) && rhs.type() == typeid(Reals)) {
    promote_matrix<dim>(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Matrix<dim, dim>)) {
    promote_matrix<dim>(size, rhs);
  }
}

template <Int dim>
void promote(LO size, Teuchos::any& lhs, Teuchos::any& rhs) {
  if (lhs.type() == rhs.type()) return;
  promote_bools(size, lhs, rhs);
  promote_scalars(size, lhs, rhs);
  promote_vectors<dim>(size, lhs, rhs);
  promote_matrices<dim>(size, lhs, rhs);
}

void promote(LO size, Int dim, Teuchos::any& lhs, Teuchos::any& rhs) {
  if (dim == 3) promote<3>(size, lhs, rhs);
  else if (dim == 2) promote<2>(size, lhs, rhs);
  else if (dim == 1) promote<1>(size, lhs, rhs);
  OMEGA_H_NORETURN();
}

template <Int dim>
void promote(LO size, Teuchos::any& x) {
  promote_bool(size, x);
  promote_scalar(size, x);
  promote_vector<dim>(size, x);
  promote_matrix<dim>(size, x);
}

void promote(LO size, Int dim, Teuchos::any& x) {
  if (dim == 3) promote<3>(size, x);
  else if (dim == 2) promote<2>(size, x);
  else if (dim == 1) promote<1>(size, x);
  OMEGA_H_NORETURN();
}

template <Int dim>
Teuchos::any add(Teuchos::any& lhs, Teuchos::any& rhs) {
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

Teuchos::any add(Int dim, Teuchos::any& lhs, Teuchos::any& rhs) {
  if (dim == 3) return add<3>(lhs, rhs);
  if (dim == 2) return add<2>(lhs, rhs);
  if (dim == 1) return add<1>(lhs, rhs);
  OMEGA_H_NORETURN(Teuchos::any);
}

template <Int dim>
Teuchos::any sub(Teuchos::any& lhs, Teuchos::any& rhs) {
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

Teuchos::any sub(Int dim, Teuchos::any& lhs, Teuchos::any& rhs) {
  if (dim == 3) return sub<3>(lhs, rhs);
  if (dim == 2) return sub<2>(lhs, rhs);
  if (dim == 1) return sub<1>(lhs, rhs);
  OMEGA_H_NORETURN(Teuchos::any);
}

template <Int dim>
Teuchos::any mul(Teuchos::any& lhs, Teuchos::any& rhs) {
  if (lhs.type() == typeid(Real) && rhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) * any_cast<Real>(rhs);
  /* begin multiply non-scalar by scalar (commutative) */
  } else if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Real)) {
    result = any_cast<Vector<dim>>(lhs) * any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Real) && rhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Real>(lhs) * any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Matrix<dim,dim>) && rhs.type() == typeid(Real)) {
    result = any_cast<Matrix<dim,dim>>(lhs) * any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Real) && rhs.type() == typeid(Matrix<dim,dim>)) {
    result = any_cast<Real>(lhs) * any_cast<Matrix<dim,dim>>(rhs);
  /* dot product */
  } else if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Vector<dim>>(lhs) * any_cast<Vector<dim>>(rhs);
  /* matrix * vector (non-commutative) */
  } else if (lhs.type() == typeid(Matrix<dim,dim>) && rhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Matrix<dim,dim>>(lhs) * any_cast<Vector<dim>>(rhs);
  } else {
    throw Teuchos::ParserFail("Invalid operand types to * operator");
  }
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

void ExprReader::at_shift(Teuchos::any& result_any, int token, std::string& text) override final {
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

void ExprReader::at_reduce(Teuchos::any& result, int token, std::string& text) override final {
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
      swap(result, rhs.at(0));
      break;
    case Teuchos::MathExpr::PROD_TERNARY:
      promote(size, dim, rhs.at(3), rhs.at(6));
      if (rhs.at(0).type() == typeid(bool)) {
        if (rhs.at(3).type() == typeid(Real)) {
          result = any_cast<bool>(rhs.at(0)) ?
                   any_cast<Real>(rhs.at(3)) :
                   any_cast<Real>(rhs.at(6));
        } else if (rhs.at(3).type() == typeid(Reals)) {
          result = any_cast<bool>(rhs.at(0)) ?
                   any_cast<Reals>(rhs.at(3)) :
                   any_cast<Reals>(rhs.at(6));
        } else {
          throw Teuchos::ParserFail("Invalid true value type in ternary operator");
        }
      } else if (rhs.at(0).type() == typeid(Bytes)) {
        promote(size, dim, rhs.at(3));
        promote(size, dim, rhs.at(6));
        result = Reals(ternary_each(any_cast<Bytes>(rhs.at(0)),
                   any_cast<Reals>(rhs.at(3)),
                   any_cast<Reals>(rhs.at(6))));
      } else {
        throw Teuchos::ParserFail("Invalid condition value type in ternary operator");
      }
      break;
    case Teuchos::MathExpr::PROD_OR:
      promote(size, dim, rhs.at(0), rhs.at(3));
      if (rhs.at(0).type() == typeid(bool)) {
        result = any_cast<bool>(rhs.at(0)) || any_cast<bool>(rhs.at(3));
      } else if (rhs.at(0).type() == typeid(Bytes)) {
        result = lor_each(any_cast<Bytes>(rhs.at(0)), any_cast<Bytes>(rhs.at(3)));
      } else {
        throw Teuchos::ParserFail("Invalid operand types to || operator");
      }
      break;
    case Teuchos::MathExpr::PROD_AND:
      promote(size, dim, rhs.at(0), rhs.at(3));
      if (rhs.at(0).type() == typeid(bool)) {
        result = any_cast<bool>(rhs.at(0)) && any_cast<bool>(rhs.at(3));
      } else if (rhs.at(0).type() == typeid(Bytes)) {
        result = land_each(any_cast<Bytes>(rhs.at(0)), any_cast<Bytes>(rhs.at(3)));
      } else {
        throw Teuchos::ParserFail("Invalid operand types to && operator");
      }
      break;
    case Teuchos::MathExpr::PROD_GT:
      promote(size, dim, rhs.at(0), rhs.at(3));
      if (rhs.at(0).type() == typeid(Real)) {
        result = any_cast<Real>(rhs.at(0)) > any_cast<Real>(rhs.at(3));
      } else if (rhs.at(0).type() == typeid(Reals)) {
        result = gt_each(any_cast<Bytes>(rhs.at(0)), any_cast<Bytes>(rhs.at(3)));
      } else {
        throw Teuchos::ParserFail("Invalid operand types to > operator");
      }
      break;
    case Teuchos::MathExpr::PROD_LT:
      promote(size, dim, rhs.at(0), rhs.at(3));
      if (rhs.at(0).type() == typeid(Real)) {
        result = any_cast<Real>(rhs.at(0)) < any_cast<Real>(rhs.at(3));
      } else if (rhs.at(0).type() == typeid(Reals)) {
        result = lt_each(any_cast<Bytes>(rhs.at(0)), any_cast<Bytes>(rhs.at(3)));
      } else {
        throw Teuchos::ParserFail("Invalid operand types to > operator");
      }
      break;
    case Teuchos::MathExpr::PROD_GEQ:
    case Teuchos::MathExpr::PROD_LEQ:
      throw Teuchos::ParserFail("Operators <= and >= not supported yet");
      break;
    case Teuchos::MathExpr::PROD_EQ:
      promote(size, dim, rhs.at(0), rhs.at(3));
      if (rhs.at(0).type() == typeid(Real)) {
        result = any_cast<Real>(rhs.at(0)) == any_cast<Real>(rhs.at(3));
      } else if (rhs.at(0).type() == typeid(Reals)) {
        result = eq_each(any_cast<Bytes>(rhs.at(0)), any_cast<Bytes>(rhs.at(3)));
      } else {
        throw Teuchos::ParserFail("Invalid operand types to == operator");
      }
      break;
    case Teuchos::MathExpr::PROD_ADD:
      promote(size, dim, rhs.at(0), rhs.at(3));
      add(dim, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_SUB:
      promote(size, dim, rhs.at(0), rhs.at(3));
      sub(dim, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_MUL:
      result = any_cast<double>(rhs.at(0)) * any_cast<double>(rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_DIV:
      result = any_cast<double>(rhs.at(0)) / any_cast<double>(rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_POW:
      result = std::pow(any_cast<double>(rhs.at(0)), any_cast<double>(rhs.at(3)));
      break;
    case Teuchos::MathExpr::PROD_CALL: {
      std::string& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      CallArgs& args = Teuchos::any_ref_cast<CallArgs>(rhs.at(4));
      TEUCHOS_TEST_FOR_EXCEPTION(args.n < 1 || args.n > 2, Teuchos::ParserFail,
          "Only unary and binary functions supported!\n");
      if (args.n == 1) {
        TEUCHOS_TEST_FOR_EXCEPTION(!unary_map.count(name), Teuchos::ParserFail,
            "Unknown unary function name \"" << name << "\"\n");
        Unary fptr = unary_map[name];
        result = (*fptr)(args.a0);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(!binary_map.count(name), Teuchos::ParserFail,
            "Unknown binary function name \"" << name << "\"\n");
        Binary fptr = binary_map[name];
        result = (*fptr)(args.a0, args.a1);
      }
      break;
    }
    case Teuchos::MathExpr::PROD_NO_ARGS: {
      CallArgs& args = Teuchos::make_any_ref<CallArgs>(result);
      args.n = 0;
      break;
    }
    case Teuchos::MathExpr::PROD_FIRST_ARG: {
      CallArgs& args = Teuchos::make_any_ref<CallArgs>(result);
      args.a0 = any_cast<double>(rhs.at(0));
      args.n = 1;
      break;
    }
    case Teuchos::MathExpr::PROD_NEXT_ARG: {
      CallArgs& args = Teuchos::any_ref_cast<CallArgs>(rhs.at(0));
      args.a1 = any_cast<double>(rhs.at(3));
      args.n = 2;
      swap(result, rhs.at(0));
      break;
    }
    case Teuchos::MathExpr::PROD_NEG:
      result = - any_cast<double>(rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_PARENS:
      result = any_cast<double>(rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_CONST:
      result = any_cast<double>(rhs.at(0));
      break;
    case Teuchos::MathExpr::PROD_VAR:
      throw Teuchos::ParserFail("Variables not supported!\n");
      break;
  }
}

}  // end namespace Omega_h
