Reals solve_laplacian(Mesh* mesh, Reals initial, Int width, Real tol) {
  CHECK(mesh.owners_have_all_upward(VERT));
  CHECK(initial.size() == mesh.nverts() * width);
  auto comm = mesh.comm();
  auto state = initial;
  auto star = mesh.ask_star(VERT);
  auto interior = mark_by_class_dim(mesh, VERT, mesh.dim());
  auto boundary = invert_marks(interior);
  auto b2v = collect_marked(boundary);
  auto weights = Reals(star.ab2b.size(), 1.0);
  auto bc_data = unmap(b2v, initial, width);
  Real maxdiff;
  Int niters = 0;
  do {
    auto new_state_nobc = graph_weighted_average(star, weights, state, width);
    auto new_state_w = deep_copy(new_state_nobc);
    map_into(bc_data, b2v, new_state_w, width);
    auto new_state = Reals(new_state_w);
    new_state = mesh.sync_array(VERT, new_state, width);
    auto diff = subtract_each(new_state, state);
    auto diffsq = multiply_each(diff, diff);
    auto maxdiffsq = comm->allreduce(max(diffsq), MAX);
    maxdiff = sqrt(maxdiffsq);
    state = new_state;
    ++niters;
  } while (maxdiff > tol);
  std::cout << "laplacian solve took " << niters << " iterations\n";
  return state;
}
