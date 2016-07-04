static void goal_stats(Mesh* mesh, char const* name, Int ent_dim, Reals values,
                       Real floor, Real ceil, Real minval, Real maxval) {
  auto low_marks = each_lt(values, floor);
  auto high_marks = each_gt(values, ceil);
  auto nlow = count_owned_marks(mesh, ent_dim, low_marks);
  auto nhigh = count_owned_marks(mesh, ent_dim, high_marks);
  auto ntotal = mesh->nglobal_ents(ent_dim);
  auto nmid = ntotal - nlow - nhigh;
  if (mesh->comm()->rank() == 0) {
    auto precision_before = std::cout.precision();
    std::ios::fmtflags stream_state(std::cout.flags());
    std::cout << std::fixed << std::setprecision(2);
    std::cout << ntotal << " " << plural_names[ent_dim];
    std::cout << ", " << name << " [" << minval << "," << maxval << "]";
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
    std::cout.precision(precision_before);
  }
}

static void get_minmax(Mesh* mesh, Reals values, Real* p_minval,
                       Real* p_maxval) {
  *p_minval = mesh->comm()->allreduce(min(values), OSH_MIN);
  *p_maxval = mesh->comm()->allreduce(max(values), OSH_MAX);
}

static void adapt_summary(Mesh* mesh, Real qual_floor, Real qual_ceil,
                          Real len_floor, Real len_ceil, Real minqual,
                          Real maxqual, Real minlen, Real maxlen) {
  goal_stats(mesh, "quality", mesh->dim(), mesh->ask_qualities(), qual_floor,
             qual_ceil, minqual, maxqual);
  goal_stats(mesh, "length", EDGE, mesh->ask_lengths(), len_floor, len_ceil,
             minlen, maxlen);
}

bool adapt_check(Mesh* mesh, Real qual_floor, Real qual_ceil, Real len_floor,
                 Real len_ceil, bool verbose) {
  Real minqual, maxqual;
  get_minmax(mesh, mesh->ask_qualities(), &minqual, &maxqual);
  Real minlen, maxlen;
  get_minmax(mesh, mesh->ask_lengths(), &minlen, &maxlen);
  if (qual_ceil <= minqual && len_floor <= minlen && maxlen <= len_ceil) {
    if (verbose && mesh->comm()->rank() == 0) {
      std::cout << "mesh is good: quality [" << minqual << "," << maxqual
                << "], length [" << minlen << "," << maxlen << "]\n";
    }
    return true;
  }
  if (verbose) {
    adapt_summary(mesh, qual_floor, qual_ceil, len_floor, len_ceil, minqual,
                  maxqual, minlen, maxlen);
  }
  return false;
}

bool adapt(Mesh* mesh, Real qual_floor, Real qual_ceil, Real len_floor,
           Real len_ceil, Int nlayers, Int verbosity) {
  Now t0 = now();
  auto comm = mesh->comm();
  CHECK(0.0 <= qual_floor);
  CHECK(qual_floor <= qual_ceil);
  CHECK(qual_ceil <= 1.0);
  CHECK(0.0 < len_floor);
  CHECK(2.0 * len_floor <= len_ceil);
  if (verbosity >= 1 && comm->rank() == 0) {
    std::cout << "before adapting:\n";
  }
  if (adapt_check(mesh, qual_floor, qual_ceil, len_floor, len_ceil,
                  (verbosity >= 1))) {
    return false;
  }
  auto input_qual = mesh->min_quality();
  CHECK(input_qual > 0.0);
  auto allow_qual = min2(qual_floor, input_qual);
  Now t1 = now();
  if ((verbosity >= 2) && comm->rank() == 0) {
    std::cout << "addressing edge lengths\n";
  }
  while (refine_by_size(mesh, len_ceil, allow_qual, (verbosity >= 2))) {
    if (verbosity >= 2) {
      adapt_check(mesh, qual_floor, qual_ceil, len_floor, len_ceil);
    }
  }
  while (coarsen_by_size(mesh, len_floor, allow_qual, (verbosity >= 2))) {
    if (verbosity >= 2) {
      adapt_check(mesh, qual_floor, qual_ceil, len_floor, len_ceil);
    }
  }
  Now t2 = now();
  bool first = true;
  while (mesh->min_quality() < qual_ceil) {
    if ((verbosity >= 2) && first && comm->rank() == 0) {
      std::cout << "addressing element qualities\n";
    }
    if (first) first = false;
    if (swap_edges(mesh, qual_ceil, nlayers, (verbosity >= 2))) {
      if (verbosity >= 2) {
        adapt_check(mesh, qual_floor, qual_ceil, len_floor, len_ceil);
      }
      continue;
    }
    if (coarsen_slivers(mesh, qual_ceil, nlayers, (verbosity >= 2))) {
      if (verbosity >= 2) {
        adapt_check(mesh, qual_floor, qual_ceil, len_floor, len_ceil);
      }
      continue;
    }
    if ((verbosity >= 1) && comm->rank() == 0) {
      std::cout << "adapt() could not satisfy quality\n";
    }
    break;
  }
  if (verbosity == 1) {
    if (comm->rank() == 0) {
      std::cout << "after adapting:\n";
    }
    adapt_check(mesh, qual_floor, qual_ceil, len_floor, len_ceil);
  }
  Now t3 = now();
  if (verbosity >= 1 && comm->rank() == 0) {
    std::cout << "addressing edge lengths took " << (t2 - t1) << " seconds\n";
  }
  if (!first && verbosity >= 1 && comm->rank() == 0) {
    std::cout << "addressing element qualities took " << (t3 - t2)
              << " seconds\n";
  }
  if (verbosity >= 1 && comm->rank() == 0) {
    std::cout << "adapting took " << (t3 - t0) << " seconds\n\n";
  }
  return true;
}
