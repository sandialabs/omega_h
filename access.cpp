template <Int dim>
Reals repeat_symm(LO n, Matrix<dim,dim> symm) {
  Write<Real> symms(n * symm_dofs(dim));
  auto f = LAMBDA(LO i) {
    set_symm(symms, i, symm);
  };
  parallel_for(n, f);
  return symms;
}

template Reals repeat_symm<3>(LO n, Matrix<3,3> symm);
template Reals repeat_symm<2>(LO n, Matrix<2,2> symm);
