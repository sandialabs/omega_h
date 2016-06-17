template <typename T>
struct CompareArrays {
  static bool compare(Read<T> a, Read<T> b, Real, Real) {
    return a == b;
  }
};

template <>
struct CompareArrays<Real> {
  static bool compare(Read<Real> a, Read<Real> b, Real tol, Real floor) {
    return are_close(a, b, tol, floor);
  }
};

template <typename T>
static bool compare_copy_data(
    Read<T> a_data, Dist a_dist,
    Read<T> b_data, Dist b_dist,
    Int ncomps, Real tol, Real floor) {
  auto a_lin_data = reduce_data_to_owners(a_data, a_dist, ncomps);
  auto b_lin_data = reduce_data_to_owners(b_data, b_dist, ncomps);
  CHECK(a_lin_data.size() == b_lin_data.size());
  bool local_result = CompareArrays<T>::compare(
      a_lin_data, b_lin_data, tol, floor);
  auto comm = a_dist.parent_comm();
  return comm->reduce_and(local_result);
}

static Read<GO> get_local_conn(Mesh* mesh, Int dim) {
  auto h2l = mesh->ask_down(dim, dim - 1);
  auto l_globals = mesh->ask_globals(dim - 1);
  auto hl2l_globals = unmap(h2l.ab2b, l_globals, 1);
  return hl2l_globals;
}

bool compare_meshes(Mesh* a, Mesh* b, Real tol, Real floor,
    bool accept_superset, bool verbose) {
  CHECK(a->comm()->size() == b->comm()->size());
  CHECK(a->comm()->rank() == b->comm()->rank());
  auto comm = a->comm();
  auto should_print = verbose && (comm->rank() == 0);
  if (a->dim() != b->dim()) {
    if (should_print) std::cout << "mesh dimensions differ\n";
    return false;
  }
  for (Int dim = 0; dim <= a->dim(); ++dim) {
    if (a->nglobal_ents(dim) != b->nglobal_ents(dim)) {
      if (should_print) {
        std::cout << "global " << singular_names[dim] << " counts differ\n";
      }
      return false;
    }
    auto a_globals = a->ask_globals(dim);
    auto b_globals = b->ask_globals(dim);
    auto a_dist = copies_to_linear_owners(comm, a_globals);
    auto b_dist = copies_to_linear_owners(comm, b_globals);
    if (dim > 0) {
      auto a_conn = get_local_conn(a, dim);
      auto b_conn = get_local_conn(b, dim);
      auto ok = compare_copy_data(
          a_conn, a_dist,
          b_conn, b_dist,
          dim + 1, 0.0, 0.0);
      if (!ok) {
        if (should_print) {
          std::cout << singular_names[dim] << " connectivity doesn't match\n";
        }
        return false;
      }
    }
    for (Int i = 0; i < a->ntags(dim); ++i) {
      auto tag = a->get_tag(dim, i);
      auto name = tag->name();
      if (!b->has_tag(dim, name)) {
        if (should_print) {
          std::cout << singular_names[dim] << " tag \"" << name
                    << "\" exists in first mesh but not second\n";
        }
        return false;
      }
      auto ncomps = tag->ncomps();
      bool ok;
      switch (tag->type()) {
      case OSH_I8:
        ok = compare_copy_data(
            a->get_array<I8>(dim, name), a_dist,
            b->get_array<I8>(dim, name), b_dist,
            ncomps, tol, floor);
        break;
      case OSH_I32:
        ok = compare_copy_data(
            a->get_array<I32>(dim, name), a_dist,
            b->get_array<I32>(dim, name), b_dist,
            ncomps, tol, floor);
        break;
      case OSH_I64:
        ok = compare_copy_data(
            a->get_array<I64>(dim, name), a_dist,
            b->get_array<I64>(dim, name), b_dist,
            ncomps, tol, floor);
        break;
      case OSH_F64:
        ok = compare_copy_data(
            a->get_array<Real>(dim, name), a_dist,
            b->get_array<Real>(dim, name), b_dist,
            ncomps, tol, floor);
        break;
      }
      if (!ok) {
        if (should_print) {
          std::cout << singular_names[dim] << " tag \""
            << name << "\" values are different\n";
        }
        return false;
      }
    }
    if (!accept_superset) {
      for (Int i = 0; i < b->ntags(dim); ++i) {
        auto tag = b->get_tag(dim, i);
        if (!a->has_tag(dim, tag->name())) {
          if (should_print) {
            std::cout << singular_names[dim] << " tag \"" << tag->name()
              << "\" exists in second mesh but not in first\n";
          }
          return false;
        }
      }
    }
  }
  return true;
}
