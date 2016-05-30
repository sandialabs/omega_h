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

Reals vectors_2d_to_3d(Reals vecs2) {
  CHECK(vecs2.size() % 2 == 0);
  Int np = vecs2.size() / 2;
  Write<Real> vecs3(np * 3);
  auto f = LAMBDA(Int i) {
    vecs3[i * 3 + 0] = vecs2[i * 2 + 0];
    vecs3[i * 3 + 1] = vecs2[i * 2 + 1];
    vecs3[i * 3 + 2] = 0.0;
  };
  parallel_for(np, f);
  return vecs3;
}
