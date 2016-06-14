static void goal_stats(Mesh* mesh,
    char const* name,
    Int ent_dim,
    Reals values,
    Real ceil,
    Real floor,
    Real minval,
    Real maxval) {
  auto low_marks = each_lt(values, floor);
  auto high_marks = each_gt(values, ceil);
  if (mesh->could_be_shared(mesh->dim())) {
    auto own_marks = mesh->owned(mesh->dim());
    low_marks = land_each(low_marks, own_marks);
    high_marks = land_each(high_marks, own_marks);
  }
  auto nlow = mesh->comm()->allreduce(GO(sum(low_marks)), SUM);
  auto nhigh = mesh->comm()->allreduce(GO(sum(high_marks)), SUM);
  auto ntotal = mesh->nglobal_ents(mesh->dim());
  auto nmid = ntotal - nlow - nhigh;
  if (mesh->comm()->rank() == 0) {
    std::ios::fmtflags stream_state(std::cout.flags());
    std::cout << std::fixed << std::setprecision(3);
    std::cout << ntotal << " " << plural_names[ent_dim];
    std::cout << ", " << name << " " << minval << "~" << maxval;
    if (nlow) {
      std::cout << ", " << nlow << " <" << floor;
    }
    if (nmid) {
      std::cout << ", " << nmid << " in [" << floor << "," << ceil << "]";
    }
    if (nhigh) {
      std::cout << ", " << nhigh << " >" << ceil;
    }
    std::cout << '\n';
    std::cout.flags(stream_state);
  }
}

static void get_minmax(
    Mesh* mesh,
    Reals values,
    Real* p_minval,
    Real* p_maxval) {
  *p_minval = mesh->comm()->allreduce(min(values), MIN);
  *p_maxval = mesh->comm()->allreduce(max(values), MAX);
}

static void adapt_summary(Mesh* mesh,
    Real qual_floor,
    Real qual_ceil,
    Real len_floor,
    Real len_ceil,
    Real minqual,
    Real maxqual,
    Real minlen,
    Real maxlen) {
  goal_stats(mesh, "quality", mesh->dim(), mesh->ask_qualities(),
      qual_floor, qual_ceil, minqual, maxqual);
  goal_stats(mesh, "length", EDGE, mesh->ask_edge_lengths(),
      len_floor, len_ceil, minlen, maxlen);
}

bool adapt_check(Mesh* mesh,
    Real qual_floor,
    Real qual_ceil,
    Real len_floor,
    Real len_ceil) {
  Real minqual, maxqual;
  get_minmax(mesh, mesh->ask_qualities(), &minqual, &maxqual);
  Real minlen, maxlen;
  get_minmax(mesh, mesh->ask_edge_lengths(), &minlen, &maxlen);
  if (qual_ceil <= minqual &&
      len_floor <= minlen &&
      maxlen <= len_ceil) {
    if (mesh->comm()->rank() == 0) {
      std::cout << "mesh is good: quality " << minqual << "~" << maxqual
        << ", length " << minlen << "~" << maxlen << '\n';
    }
    return true;
  }
  adapt_summary(mesh,
      qual_floor, qual_ceil, len_floor, len_ceil,
      minqual, maxqual, minlen, maxlen);
  return false;
}

bool adapt(Mesh* mesh,
    Real qual_ceil,
    Real qual_floor,
    Real len_floor,
    Real len_ceil) {
  if (adapt_check(mesh, qual_ceil, qual_floor, len_floor, len_ceil)) {
    return false;
  }
  /* TODO: cascade into modification */
  return false;
}
