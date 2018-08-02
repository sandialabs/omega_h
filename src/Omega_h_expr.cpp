#include <Omega_h_expr.hpp>
#include <Omega_h_matrix.hpp>

#include <Omega_h_array_ops.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_vector.hpp>

namespace Omega_h {

using Teuchos::any;
using Teuchos::any_cast;

namespace {

void promote_bool(LO size, any& x) { x = Bytes(size, Byte(any_cast<bool>(x))); }

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

void promote_scalar(LO size, any& x) { x = Reals(size, any_cast<Real>(x)); }

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
  x = repeat_matrix(size, any_cast<Matrix<dim, dim>>(x));
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
  ScopedTimer("promote(size, dim, lhs, rhs)");
  if (dim == 3)
    promote<3>(size, lhs, rhs);
  else if (dim == 2)
    promote<2>(size, lhs, rhs);
  else if (dim == 1)
    promote<1>(size, lhs, rhs);
  else {
    OMEGA_H_NORETURN();
  }
}

template <Int dim>
void promote(LO size, any& x) {
  if (x.type() == typeid(bool)) {
    promote_bool(size, x);
  } else if (x.type() == typeid(Real)) {
    promote_scalar(size, x);
  } else if (x.type() == typeid(Vector<dim>)) {
    promote_vector<dim>(size, x);
  } else if (x.type() == typeid(Matrix<dim, dim>)) {
    promote_matrix<dim>(size, x);
  } else if (x.type() == typeid(Reals)) {
    return;
  } else {
    std::string msg;
    msg += "Unexpected type ";
    msg += x.type().name();
    msg += " being promoted";
    throw Teuchos::ParserFail(msg);
  }
}

void promote(LO size, Int dim, any& x) {
  ScopedTimer("promote(size, dim, any)");
  if (dim == 3)
    promote<3>(size, x);
  else if (dim == 2)
    promote<2>(size, x);
  else if (dim == 1)
    promote<1>(size, x);
  else
    OMEGA_H_NORETURN();
}

void ternary(
    LO size, Int dim, any& result, any& cond, any& true_val, any& false_val) {
  OMEGA_H_TIME_FUNCTION;
  if (cond.type() == typeid(bool)) {
    if (true_val.type() == typeid(Real)) {
      result = any_cast<bool>(cond) ? any_cast<Real>(true_val)
                                    : any_cast<Real>(false_val);
    } else if (true_val.type() == typeid(Reals)) {
      result = any_cast<bool>(cond) ? any_cast<Reals>(true_val)
                                    : any_cast<Reals>(false_val);
    } else {
      throw Teuchos::ParserFail("Invalid true value type in ternary operator");
    }
  } else if (cond.type() == typeid(Bytes)) {
    promote(size, dim, true_val);
    promote(size, dim, false_val);
    result = Reals(ternary_each(any_cast<Bytes>(cond),
        any_cast<Reals>(true_val), any_cast<Reals>(false_val)));
  } else {
    throw Teuchos::ParserFail(
        "Invalid condition value type in ternary operator");
  }
}

void eval_or(any& result, any& lhs, any& rhs) {
  if (lhs.type() == typeid(bool)) {
    result = any_cast<bool>(lhs) || any_cast<bool>(rhs);
  } else if (lhs.type() == typeid(Bytes)) {
    result = lor_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to || operator");
  }
}

void eval_and(any& result, any& lhs, any& rhs) {
  if (lhs.type() == typeid(bool)) {
    result = any_cast<bool>(lhs) && any_cast<bool>(rhs);
  } else if (lhs.type() == typeid(Bytes)) {
    result = land_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to && operator");
  }
}

void gt(any& result, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) > any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = gt_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to > operator");
  }
}

void lt(any& result, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) < any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = lt_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to > operator");
  }
}

void eq(any& result, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) == any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = eq_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to == operator");
  }
}

template <Int dim>
void add(any& result, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) + any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Vector<dim>>(lhs) + any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Matrix<dim, dim>)) {
    result = any_cast<Matrix<dim, dim>>(lhs) + any_cast<Matrix<dim, dim>>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = Reals(add_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to + operator");
  }
}

void add(Int dim, any& result, any& lhs, any& rhs) {
  if (dim == 3) return add<3>(result, lhs, rhs);
  if (dim == 2) return add<2>(result, lhs, rhs);
  if (dim == 1) return add<1>(result, lhs, rhs);
}

template <Int dim>
void sub(any& result, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) - any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Vector<dim>>(lhs) - any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Matrix<dim, dim>)) {
    result = any_cast<Matrix<dim, dim>>(lhs) - any_cast<Matrix<dim, dim>>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    result = Reals(subtract_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else {
    throw Teuchos::ParserFail("Invalid operand types to + operator");
  }
}

void sub(Int dim, any& result, any& lhs, any& rhs) {
  if (dim == 3) return sub<3>(result, lhs, rhs);
  if (dim == 2) return sub<2>(result, lhs, rhs);
  if (dim == 1) return sub<1>(result, lhs, rhs);
}

template <Int dim>
void mul(LO size, any& result, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real) && rhs.type() == typeid(Real)) {
    result = any_cast<Real>(lhs) * any_cast<Real>(rhs);
    /* begin multiply non-scalar by scalar (commutative) */
  } else if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Real)) {
    result = any_cast<Vector<dim>>(lhs) * any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Real) && rhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Real>(lhs) * any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Matrix<dim, dim>) &&
             rhs.type() == typeid(Real)) {
    result = any_cast<Matrix<dim, dim>>(lhs) * any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Real) &&
             rhs.type() == typeid(Matrix<dim, dim>)) {
    result = any_cast<Real>(lhs) * any_cast<Matrix<dim, dim>>(rhs);
    /* dot product */
  } else if (lhs.type() == typeid(Vector<dim>) &&
             rhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Vector<dim>>(lhs) * any_cast<Vector<dim>>(rhs);
    /* matrix * vector (non-commutative) */
  } else if (lhs.type() == typeid(Matrix<dim, dim>) &&
             rhs.type() == typeid(Vector<dim>)) {
    result = any_cast<Matrix<dim, dim>>(lhs) * any_cast<Vector<dim>>(rhs);
    /* matrix * matrix (non-commutative) */
  } else if (lhs.type() == typeid(Matrix<dim, dim>) &&
             rhs.type() == typeid(Matrix<dim, dim>)) {
    result = any_cast<Matrix<dim, dim>>(lhs) * any_cast<Matrix<dim, dim>>(rhs);
  } else if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Reals)) {
    auto& lhs_vals = any_cast<Reals>(lhs);
    auto& rhs_vals = any_cast<Reals>(rhs);
    if (rhs_vals.size() == size) {  // RHS is scalars
      result = Reals(multiply_each(lhs_vals, rhs_vals));
    } else if (lhs_vals.size() == size) {  // LHS is scalars
      result = Reals(multiply_each(rhs_vals, lhs_vals));
    } else if (lhs_vals.size() == size * dim &&
               rhs_vals.size() == size * dim) {  // dot products
      result = dot_vectors(lhs_vals, rhs_vals, dim);
    } else if (lhs_vals.size() == size * matrix_ncomps(dim, dim) &&
               rhs_vals.size() == size * dim) {  // matrices * vectors
      result = matrices_times_vectors(lhs_vals, rhs_vals, dim);
    } else if (lhs_vals.size() == size * matrix_ncomps(dim, dim) &&
               rhs_vals.size() == size * matrix_ncomps(dim, dim)) {
      result = matrices_times_matrices(lhs_vals, rhs_vals, dim);
    } else {
      throw Teuchos::ParserFail("Unexpected array size in * operator");
    }
  } else {
    throw Teuchos::ParserFail("Invalid operand types to * operator");
  }
}

void mul(LO size, Int dim, any& result, any& lhs, any& rhs) {
  if (dim == 3) mul<3>(size, result, lhs, rhs);
  if (dim == 2) mul<2>(size, result, lhs, rhs);
  if (dim == 1) mul<1>(size, result, lhs, rhs);
}

template <Int dim>
void div(any& result, any& lhs, any& rhs) {
  if (rhs.type() == typeid(Reals)) {
    result = Reals(divide_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else if (rhs.type() == typeid(Real)) {
    if (lhs.type() == typeid(Real)) {
      result = any_cast<Real>(lhs) / any_cast<Real>(rhs);
    } else if (lhs.type() == typeid(Vector<dim>)) {
      result = any_cast<Vector<dim>>(lhs) / any_cast<Real>(rhs);
    } else if (lhs.type() == typeid(Matrix<dim, dim>)) {
      result = any_cast<Matrix<dim, dim>>(lhs) / any_cast<Real>(rhs);
    } else {
      throw Teuchos::ParserFail("Invalid left operand type in / operator");
    }
  } else {
    throw Teuchos::ParserFail("Invalid right operand type in / operator");
  }
}

void div(Int dim, any& result, any& lhs, any& rhs) {
  if (dim == 3) div<3>(result, lhs, rhs);
  if (dim == 2) div<2>(result, lhs, rhs);
  if (dim == 1) div<1>(result, lhs, rhs);
}

template <Int dim>
void eval_pow(any& result, any& lhs, any& rhs) {
  if (rhs.type() == typeid(Reals)) {
    result = Reals(pow_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else if (rhs.type() == typeid(Real)) {
    result = std::pow(any_cast<Real>(lhs), any_cast<Real>(rhs));
  } else {
    throw Teuchos::ParserFail("Invalid right operand type in ^ operator");
  }
}

void eval_pow(Int dim, any& result, any& lhs, any& rhs) {
  if (dim == 3) eval_pow<3>(result, lhs, rhs);
  if (dim == 2) eval_pow<2>(result, lhs, rhs);
  if (dim == 1) eval_pow<1>(result, lhs, rhs);
}

template <Int dim>
void neg(any& result, any& val) {
  if (val.type() == typeid(Real)) {
    result = -any_cast<Real>(val);
  } else if (val.type() == typeid(Vector<dim>)) {
    result = -any_cast<Vector<dim>>(val);
  } else if (val.type() == typeid(Matrix<dim, dim>)) {
    result = -any_cast<Matrix<dim, dim>>(val);
  } else if (val.type() == typeid(Reals)) {
    result = Reals(multiply_each_by(any_cast<Reals>(val), -1.0));
  } else {
    throw Teuchos::ParserFail("Invalid operand type to negation operator");
  }
}

void neg(Int dim, any& result, any& val) {
  if (dim == 3) neg<3>(result, val);
  if (dim == 2) neg<2>(result, val);
  if (dim == 1) neg<1>(result, val);
}

template <Int dim>
void access(LO size, any& result, any& var, ExprReader::Args& args) {
  auto i = static_cast<Int>(any_cast<Real>(args.at(0)));
  auto j =
      args.size() > 1 ? static_cast<Int>(any_cast<Real>(args.at(0))) : Int(-1);
  if (var.type() == typeid(Vector<dim>)) {
    result = (any_cast<Vector<dim>>(var))(i);
  } else if (var.type() == typeid(Matrix<dim, dim>)) {
    result = (any_cast<Matrix<dim, dim>>(var))(i, j);
  } else if (var.type() == typeid(Reals)) {
    auto array = any_cast<Reals>(var);
    if (array.size() == size * dim) {
      result = Reals(get_component(array, dim, i));
    } else if (array.size() == size * matrix_ncomps(dim, dim)) {
      result =
          Reals(get_component(array, matrix_ncomps(dim, dim), j * dim + i));
    } else {
      throw Teuchos::ParserFail("Unexpected array size in access operator\n");
    }
  } else {
    throw Teuchos::ParserFail("Unexpected variable type in access operator\n");
  }
}

void access(LO size, Int dim, any& result, any& var, ExprReader::Args& args) {
  if (dim == 3) access<3>(size, result, var, args);
  if (dim == 2) access<2>(size, result, var, args);
  if (dim == 1) access<1>(size, result, var, args);
}

template <Int dim>
void make_vector(any& result, ExprReader::Args& args) {
  auto v = zero_vector<dim>();
  Int i;
  for (i = 0; i < Int(args.size()); ++i) {
    auto& arg = args[std::size_t(i)];
    v[i] = any_cast<Real>(arg);
  }
  for (; i < dim; ++i) {
    v[i] = v[Int(args.size() - 1)];
  }
  result = v;
}

void make_vector(LO size, Int dim, any& result, ExprReader::Args& args) {
  TEUCHOS_TEST_FOR_EXCEPTION(args.size() > std::size_t(dim),
      Teuchos::ParserFail, "Too many arguments to vector()\n");
  bool has_arrays = false;
  for (auto& arg : args)
    if (arg.type() == typeid(Reals)) has_arrays = true;
  if (has_arrays) {
    std::vector<Read<Real>> arrays;
    Int i;
    for (i = 0; i < Int(args.size()); ++i) {
      auto& arg = args[std::size_t(i)];
      promote(size, dim, arg);
      arrays.push_back(any_cast<Reals>(arg));
    }
    for (; i < dim; ++i) {
      arrays.push_back(arrays.back());
    }
    OMEGA_H_CHECK(Int(arrays.size()) == dim);
    result = Reals(interleave(arrays));
  } else {
    if (dim == 3) make_vector<3>(result, args);
    if (dim == 2) make_vector<2>(result, args);
    if (dim == 1) make_vector<1>(result, args);
  }
}

void eval_exp(LO size, any& result, ExprReader::Args& args) {
  TEUCHOS_TEST_FOR_EXCEPTION(args.size() != 1, Teuchos::ParserFail,
      "exp() takes exactly one argument, given " << args.size());
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    result = std::exp(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = Teuchos::any_cast<Reals>(in_any);
    TEUCHOS_TEST_FOR_EXCEPTION(a.size() != size, Teuchos::ParserFail,
        "exp() given array that wasn't scalars");
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::exp(a[i]); };
    parallel_for(a.size(), f, "eval_exp(Reals)");
    result = Reals(out);
  } else {
    throw Teuchos::ParserFail("unexpected argument type to exp()");
  }
}

void eval_sqrt(LO size, any& result, ExprReader::Args& args) {
  TEUCHOS_TEST_FOR_EXCEPTION(args.size() != 1, Teuchos::ParserFail,
      "sqrt() takes exactly one argument, given " << args.size());
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    result = std::sqrt(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = Teuchos::any_cast<Reals>(in_any);
    TEUCHOS_TEST_FOR_EXCEPTION(a.size() != size, Teuchos::ParserFail,
        "sqrt() given array that wasn't scalars");
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::sqrt(a[i]); };
    parallel_for(a.size(), f, "eval_sqrt(Reals)");
    result = Reals(out);
  } else {
    throw Teuchos::ParserFail("unexpected argument type to sqrt()");
  }
}

void eval_sin(LO size, any& result, ExprReader::Args& args) {
  TEUCHOS_TEST_FOR_EXCEPTION(args.size() != 1, Teuchos::ParserFail,
      "sin() takes exactly one argument, given " << args.size());
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    result = std::sin(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = Teuchos::any_cast<Reals>(in_any);
    TEUCHOS_TEST_FOR_EXCEPTION(a.size() != size, Teuchos::ParserFail,
        "sin() given array that wasn't scalars");
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::sin(a[i]); };
    parallel_for(a.size(), f, "eval_sin(Reals)");
    result = Reals(out);
  } else {
    throw Teuchos::ParserFail("unexpected argument type to sin()");
  }
}

void eval_cos(LO size, any& result, ExprReader::Args& args) {
  TEUCHOS_TEST_FOR_EXCEPTION(args.size() != 1, Teuchos::ParserFail,
      "cos() takes exactly one argument, given " << args.size());
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    result = std::cos(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = Teuchos::any_cast<Reals>(in_any);
    TEUCHOS_TEST_FOR_EXCEPTION(a.size() != size, Teuchos::ParserFail,
        "cos() given array that wasn't scalars");
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::cos(a[i]); };
    parallel_for(a.size(), f, "eval_cos(Reals)");
    result = Reals(out);
  } else {
    throw Teuchos::ParserFail("unexpected argument type to cos()");
  }
}

template <Int dim>
void eval_norm(LO size, any& result, ExprReader::Args& args) {
  TEUCHOS_TEST_FOR_EXCEPTION(args.size() != 1, Teuchos::ParserFail,
      "norm() takes exactly one argument, given " << args.size());
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Vector<dim>)) {
    result = Omega_h::norm(any_cast<Vector<dim>>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = Teuchos::any_cast<Reals>(in_any);
    TEUCHOS_TEST_FOR_EXCEPTION(a.size() != size * dim, Teuchos::ParserFail,
        "norm() given array that wasn't vectors");
    auto out = Write<Real>(size);
    auto f = OMEGA_H_LAMBDA(LO i) {
      out[i] = Omega_h::norm(Omega_h::get_vector<dim>(a, i));
    };
    parallel_for(size, f, "eval_norm(Reals)");
    result = Reals(out);
  } else {
    throw Teuchos::ParserFail("unexpected argument type to norm()");
  }
}

void eval_norm(Int dim, LO size, any& result, ExprReader::Args& args) {
  if (dim == 3) eval_norm<3>(size, result, args);
  if (dim == 2) eval_norm<2>(size, result, args);
  if (dim == 1) eval_norm<1>(size, result, args);
}

}  // end anonymous namespace

ExprEnv::ExprEnv(LO size_in, Int dim_in) : size(size_in), dim(dim_in) {
  auto local_size = size;
  auto local_dim = dim;
  auto vector = [=](any& result, Args& args) {
    make_vector(local_size, local_dim, result, args);
  };
  register_function("exp",
      [=](any& result, Args& args) { eval_exp(local_size, result, args); });
  register_function("sqrt",
      [=](any& result, Args& args) { eval_sqrt(local_size, result, args); });
  register_function("sin",
      [=](any& result, Args& args) { eval_sin(local_size, result, args); });
  register_function("cos",
      [=](any& result, Args& args) { eval_cos(local_size, result, args); });
  register_function("vector", vector);
  register_function("norm", [=](any& result, Args& args) {
    eval_norm(local_dim, local_size, result, args);
  });
  register_variable("d", any(Real(dim)));
  if (dim == 3) register_variable("I", any(identity_matrix<3, 3>()));
  if (dim == 2) register_variable("I", any(identity_matrix<2, 2>()));
  if (dim == 1) register_variable("I", any(identity_matrix<1, 1>()));
  register_variable("pi", any(Real(Omega_h::PI)));
}

void ExprEnv::register_variable(
    std::string const& name, Teuchos::any const& value) {
  OMEGA_H_TIME_FUNCTION;
  // OMEGA_H_CHECK(variables.find(name) == variables.end());
  OMEGA_H_CHECK(functions.find(name) == functions.end());
  variables[name] = value;
}

void ExprEnv::register_function(
    std::string const& name, Function const& value) {
  OMEGA_H_CHECK(variables.find(name) == variables.end());
  // OMEGA_H_CHECK(functions.find(name) == functions.end());
  functions[name] = value;
}

void ExprEnv::repeat(Teuchos::any& x) { promote(size, dim, x); }

ExprReader::ExprReader(LO size_in, Int dim_in)
    : Teuchos::Reader(Teuchos::MathExpr::ask_reader_tables()),
      env(size_in, dim_in) {}

ExprReader::~ExprReader() {}

void ExprReader::register_variable(
    std::string const& name, Teuchos::any const& value) {
  env.register_variable(name, value);
}

void ExprReader::register_function(
    std::string const& name, Function const& value) {
  env.register_function(name, value);
}

void ExprReader::repeat(Teuchos::any& x) { env.repeat(x); }

void ExprReader::at_shift(any& result_any, int token, std::string& text) {
  using std::swap;
  switch (token) {
    case Teuchos::MathExpr::TOK_NAME: {
      auto& result = Teuchos::make_any_ref<std::string>(result_any);
      swap(result, text);
      break;
    }
    case Teuchos::MathExpr::TOK_CONST: {
      result_any = std::atof(text.c_str());
    }
  }
}

void ExprReader::at_reduce(any& result, int prod, std::vector<any>& rhs) {
  OMEGA_H_TIME_FUNCTION;
  using std::swap;
  switch (prod) {
    case Teuchos::MathExpr::PROD_PROGRAM: {
      TEUCHOS_TEST_FOR_EXCEPTION(rhs.at(1).empty(), Teuchos::ParserFail,
          "Omega_h::ExprReader needs an expression to evaluate!");
      swap(result, rhs.at(1));
      break;
    }
    case Teuchos::MathExpr::PROD_NO_STATEMENTS:
    case Teuchos::MathExpr::PROD_NO_EXPR:
    case Teuchos::MathExpr::PROD_NEXT_STATEMENT: {
      break;
    }
    case Teuchos::MathExpr::PROD_ASSIGN: {
      std::string const& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      swap(env.variables[name], rhs.at(4));
      break;
    }
    case Teuchos::MathExpr::PROD_YES_EXPR:
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
      promote(env.size, env.dim, rhs.at(3), rhs.at(6));
      ternary(env.size, env.dim, result, rhs.at(0), rhs.at(3), rhs.at(6));
      break;
    case Teuchos::MathExpr::PROD_OR:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      eval_or(result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_AND:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      eval_and(result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_GT:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      gt(result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_LT:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      lt(result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_GEQ:
    case Teuchos::MathExpr::PROD_LEQ:
      throw Teuchos::ParserFail("Operators <= and >= not supported yet");
    case Teuchos::MathExpr::PROD_EQ:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      eq(result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_ADD:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      add(env.dim, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_SUB:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      sub(env.dim, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_MUL:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      mul(env.size, env.dim, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_DIV:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      div(env.dim, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_POW:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      eval_pow(env.dim, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_CALL: {
      ScopedTimer("call");
      auto& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      auto& args = Teuchos::any_ref_cast<Args>(rhs.at(4));
      auto vit = env.variables.find(name);
      if (vit == env.variables.end()) {
        /* function call */
        auto fit = env.functions.find(name);
        TEUCHOS_TEST_FOR_EXCEPTION(fit == env.functions.end(),
            Teuchos::ParserFail,
            "\"" << name << "\" is neither a variable nor a function name\n");
        fit->second(result, args);
      } else {
        /* access operator for vector/matrix */
        auto& val = vit->second;
        access(env.size, env.dim, result, val, args);
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
      neg(env.dim, result, rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_VAL_PARENS:
    case Teuchos::MathExpr::PROD_BOOL_PARENS:
      swap(result, rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_VAR:
      auto& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      auto it = env.variables.find(name);
      TEUCHOS_TEST_FOR_EXCEPTION(it == env.variables.end(), Teuchos::ParserFail,
          "unknown variable name \"" << name << "\"\n");
      result = it->second;
      break;
  }
}

ExprOp::~ExprOp() {}

using OpPtr = std::shared_ptr<ExprOp>;

struct ConstOp : public ExprOp {
  double value;
  OpPtr rhs;
  virtual ~ConstOp() override final = default;
  ConstOp(double value_in) : value(value_in) {}
  virtual void eval(ExprEnv& env, Teuchos::any& result) override final;
};
void ConstOp::eval(ExprEnv&, Teuchos::any& result) { result = value; }

struct AssignOp : public ExprOp {
  std::string name;
  OpPtr rhs;
  virtual ~AssignOp() override final = default;
  AssignOp(std::string const& name_in, OpPtr rhs_in)
      : name(name_in), rhs(rhs_in) {}
  virtual void eval(ExprEnv& env, Teuchos::any& result) override final;
};
void AssignOp::eval(ExprEnv& env, Teuchos::any& result) {
  using std::swap;
  rhs->eval(env, result);
  swap(env.variables[name], result);
}

struct VarOp : public ExprOp {
  std::string name;
  virtual ~VarOp() override final = default;
  VarOp(std::string const& name_in) : name(name_in) {}
  virtual void eval(ExprEnv& env, Teuchos::any& result) override final;
};
void VarOp::eval(ExprEnv& env, Teuchos::any& result) {
  auto it = env.variables.find(name);
  TEUCHOS_TEST_FOR_EXCEPTION(it == env.variables.end(), Teuchos::ParserFail,
      "unknown variable name \"" << name << "\"\n");
  result = it->second;
}

struct NegOp : public ExprOp {
  OpPtr rhs;
  virtual ~NegOp() override final = default;
  NegOp(OpPtr rhs_in) : rhs(rhs_in) {}
  virtual void eval(ExprEnv& env, Teuchos::any& result) override final;
};
void NegOp::eval(ExprEnv& env, Teuchos::any& result) {
  Teuchos::any rhs_val;
  rhs->eval(env, rhs_val);
  neg(env.dim, result, rhs_val);
}

struct TernaryOp : public ExprOp {
  OpPtr cond;
  OpPtr lhs;
  OpPtr rhs;
  virtual ~TernaryOp() override final = default;
  TernaryOp(OpPtr cond_in, OpPtr lhs_in, OpPtr rhs_in)
      : cond(cond_in), lhs(lhs_in), rhs(rhs_in) {}
  virtual void eval(ExprEnv& env, Teuchos::any& result) override final;
};
void TernaryOp::eval(ExprEnv& env, Teuchos::any& result) {
  Teuchos::any cond_val, lhs_val, rhs_val;
  cond->eval(env, cond_val);
  lhs->eval(env, lhs_val);
  rhs->eval(env, rhs_val);
  promote(env.size, env.dim, lhs_val, rhs_val);
  ternary(env.size, env.dim, result, cond_val, lhs_val, rhs_val);
}

struct CallOp : public ExprOp {
  std::string name;
  std::vector<OpPtr> rhs;
  ExprEnv::Args args;
  virtual ~CallOp() override final = default;
  CallOp(std::string const& name_in, ExprEnv::Args const& args_in)
      : name(name_in) {
    for (auto& arg : args_in) {
      OpPtr op = Teuchos::any_cast<OpPtr>(arg);
      rhs.push_back(op);
    }
    args.reserve(rhs.size());
  }
  virtual void eval(ExprEnv& env, Teuchos::any& result) override final;
};
void CallOp::eval(ExprEnv& env, Teuchos::any& result) {
  args.resize(rhs.size());
  for (std::size_t i = 0; i < rhs.size(); ++i) {
    rhs[i]->eval(env, args[i]);
  }
  auto vit = env.variables.find(name);
  if (vit == env.variables.end()) {
    /* function call */
    auto fit = env.functions.find(name);
    TEUCHOS_TEST_FOR_EXCEPTION(fit == env.functions.end(), Teuchos::ParserFail,
        "\"" << name << "\" is neither a variable nor a function name\n");
    fit->second(result, args);
  } else {
    /* access operator for vector/matrix */
    auto& val = vit->second;
    access(env.size, env.dim, result, val, args);
  }
  args.clear();
}

#define OMEGA_H_BINARY_OP(ClassName, func_call)                                \
  struct ClassName : public ExprOp {                                           \
    OpPtr lhs;                                                                 \
    OpPtr rhs;                                                                 \
    virtual ~ClassName() override final = default;                             \
    ClassName(OpPtr lhs_in, OpPtr rhs_in) : lhs(lhs_in), rhs(rhs_in) {}        \
    virtual void eval(ExprEnv& env, Teuchos::any& result) override final;      \
  };                                                                           \
  void ClassName::eval(ExprEnv& env, Teuchos::any& result) {                   \
    Teuchos::any lhs_val, rhs_val;                                             \
    lhs->eval(env, lhs_val);                                                   \
    rhs->eval(env, rhs_val);                                                   \
    promote(env.size, env.dim, lhs_val, rhs_val);                              \
    func_call;                                                                 \
  }

OMEGA_H_BINARY_OP(OrOp, eval_or(result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(AndOp, eval_and(result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(GtOp, gt(result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(LtOp, lt(result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(EqOp, eq(result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(AddOp, add(env.dim, result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(SubOp, sub(env.dim, result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(MulOp, mul(env.size, env.dim, result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(DivOp, div(env.dim, result, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(PowOp, eval_pow(env.dim, result, lhs_val, rhs_val));

#undef OMEGA_H_BINARY_OP

#define OMEGA_H_BINARY_REDUCE(ClassName)                                       \
  {                                                                            \
    OpPtr lhs_op = Teuchos::any_cast<OpPtr>(rhs.at(0));                        \
    OpPtr rhs_op = Teuchos::any_cast<OpPtr>(rhs.at(3));                        \
    OpPtr result_op(new ClassName(lhs_op, rhs_op));                            \
    result = result_op;                                                        \
  }

ExprOpsReader::ExprOpsReader()
    : Teuchos::Reader(Teuchos::MathExpr::ask_reader_tables()) {}

std::shared_ptr<Omega_h::ExprOp> ExprOpsReader::read_ops(std::string const& expr) {
  Teuchos::any any_op;
  read_string(any_op, expr, expr);
  return Teuchos::any_cast<std::shared_ptr<Omega_h::ExprOp>>(any_op);
}

void ExprOpsReader::at_shift(any& result_any, int token, std::string& text) {
  using std::swap;
  switch (token) {
    case Teuchos::MathExpr::TOK_NAME: {
      auto& result = Teuchos::make_any_ref<std::string>(result_any);
      swap(result, text);
      break;
    }
    case Teuchos::MathExpr::TOK_CONST: {
      OpPtr result_op(new ConstOp(std::atof(text.c_str())));
      result_any = result_op;
    }
  }
}

void ExprOpsReader::at_reduce(any& result, int prod, std::vector<any>& rhs) {
  OMEGA_H_TIME_FUNCTION;
  using std::swap;
  switch (prod) {
    case Teuchos::MathExpr::PROD_PROGRAM: {
      TEUCHOS_TEST_FOR_EXCEPTION(rhs.at(1).empty(), Teuchos::ParserFail,
          "Omega_h::ExprReader needs an expression to evaluate!");
      swap(result, rhs.at(1));
      break;
    }
    case Teuchos::MathExpr::PROD_NO_STATEMENTS:
    case Teuchos::MathExpr::PROD_NO_EXPR:
    case Teuchos::MathExpr::PROD_NEXT_STATEMENT: {
      break;
    }
    case Teuchos::MathExpr::PROD_ASSIGN: {
      std::string const& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      OpPtr rhs_op = Teuchos::any_cast<OpPtr>(rhs.at(4));
      OpPtr op(new AssignOp(name, rhs_op));
      result = op;
      break;
    }
    case Teuchos::MathExpr::PROD_YES_EXPR:
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
    case Teuchos::MathExpr::PROD_TERNARY: {
      OpPtr cond_op = Teuchos::any_cast<OpPtr>(rhs.at(0));
      OpPtr lhs_op = Teuchos::any_cast<OpPtr>(rhs.at(3));
      OpPtr rhs_op = Teuchos::any_cast<OpPtr>(rhs.at(6));
      OpPtr result_op(new TernaryOp(cond_op, lhs_op, rhs_op));
      result = result_op;
    } break;
    case Teuchos::MathExpr::PROD_OR:
      OMEGA_H_BINARY_REDUCE(OrOp);
      break;
    case Teuchos::MathExpr::PROD_AND:
      OMEGA_H_BINARY_REDUCE(AndOp);
      break;
    case Teuchos::MathExpr::PROD_GT:
      OMEGA_H_BINARY_REDUCE(GtOp);
      break;
    case Teuchos::MathExpr::PROD_LT:
      OMEGA_H_BINARY_REDUCE(LtOp);
      break;
    case Teuchos::MathExpr::PROD_GEQ:
    case Teuchos::MathExpr::PROD_LEQ:
      throw Teuchos::ParserFail("Operators <= and >= not supported yet");
    case Teuchos::MathExpr::PROD_EQ:
      OMEGA_H_BINARY_REDUCE(EqOp);
      break;
    case Teuchos::MathExpr::PROD_ADD:
      OMEGA_H_BINARY_REDUCE(AddOp);
      break;
    case Teuchos::MathExpr::PROD_SUB:
      OMEGA_H_BINARY_REDUCE(SubOp);
      break;
    case Teuchos::MathExpr::PROD_MUL:
      OMEGA_H_BINARY_REDUCE(MulOp);
      break;
    case Teuchos::MathExpr::PROD_DIV:
      OMEGA_H_BINARY_REDUCE(DivOp);
      break;
    case Teuchos::MathExpr::PROD_POW:
      OMEGA_H_BINARY_REDUCE(PowOp);
      break;
    case Teuchos::MathExpr::PROD_CALL: {
      auto& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      auto& args = Teuchos::any_ref_cast<ExprEnv::Args>(rhs.at(4));
      OpPtr result_op(new CallOp(name, args));
      result = result_op;
      break;
    }
    case Teuchos::MathExpr::PROD_NO_ARGS: {
      result = ExprEnv::Args{};
      break;
    }
    case Teuchos::MathExpr::PROD_FIRST_ARG: {
      auto& args = Teuchos::make_any_ref<ExprEnv::Args>(result);
      args.push_back(any());
      swap(args.back(), rhs.at(0));
      break;
    }
    case Teuchos::MathExpr::PROD_NEXT_ARG: {
      auto& args = Teuchos::any_ref_cast<ExprEnv::Args>(rhs.at(0));
      args.push_back(any());
      swap(args.back(), rhs.at(3));
      swap(result, rhs.at(0));
      break;
    }
    case Teuchos::MathExpr::PROD_NEG: {
      OpPtr rhs_op = Teuchos::any_cast<OpPtr>(rhs.at(2));
      OpPtr result_op(new NegOp(rhs_op));
      result = result_op;
    } break;
    case Teuchos::MathExpr::PROD_VAL_PARENS:
    case Teuchos::MathExpr::PROD_BOOL_PARENS:
      swap(result, rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_VAR:
      auto& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      OpPtr result_op(new VarOp(name));
      result = result_op;
      break;
  }
}

#undef OMEGA_H_BINARY_REDUCE

}  // end namespace Omega_h
