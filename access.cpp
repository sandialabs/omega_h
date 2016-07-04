template <Int dim>
Reals repeat_symm(LO n, Matrix<dim, dim> symm) {
  Write<Real> symms(n * symm_dofs(dim));
  auto f = LAMBDA(LO i) { set_symm(symms, i, symm); };
  parallel_for(n, f);
  return symms;
}

template Reals repeat_symm<3>(LO n, Matrix<3, 3> symm);
template Reals repeat_symm<2>(LO n, Matrix<2, 2> symm);

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

Reals vectors_3d_to_2d(Reals vecs3) {
  CHECK(vecs3.size() % 3 == 0);
  Int np = vecs3.size() / 3;
  Write<Real> vecs2(np * 2);
  auto f = LAMBDA(Int i) {
    vecs2[i * 2 + 0] = vecs3[i * 3 + 0];
    vecs2[i * 2 + 1] = vecs3[i * 3 + 1];
  };
  parallel_for(np, f);
  return vecs2;
}

Reals average_field(Mesh* mesh, Int dim, LOs a2e, Int ncomps, Reals v2x) {
  auto ev2v = mesh->ask_verts_of(dim);
  auto degree = simplex_degrees[dim][VERT];
  CHECK(v2x.size() % ncomps == 0);
  auto na = a2e.size();
  Write<Real> out(na * ncomps);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    for (Int j = 0; j < ncomps; ++j) {
      Real comp = 0;
      for (Int k = 0; k < degree; ++k) {
        auto v = ev2v[e * degree + k];
        comp += v2x[v * ncomps + j];
      }
      comp /= degree;
      out[a * ncomps + j] = comp;
    }
  };
  parallel_for(na, f);
  return out;
}
