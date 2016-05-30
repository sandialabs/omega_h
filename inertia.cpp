namespace inertia {

namespace {

/* since the cross product is involved, we do this all
   in 3D only. */

Vector<3> get_center(CommPtr comm, Reals coords, Reals masses,
    Real total_mass) {
  auto n = masses.size();
  Write<Real> weighted_coords(n * 3);
  auto f = LAMBDA(LO i) {
    set_vec<3>(weighted_coords, i, masses[i] * get_vec<3>(coords, i));
  };
  parallel_for(n, f);
  Vector<3> result;
  repro_sum(comm, Reals(weighted_coords), 3, &result[0]);
  return result / total_mass;
}

Matrix<3,3> get_matrix(CommPtr comm, Reals coords,
    Reals masses, Vector<3> center) {
  auto n = masses.size();
  Write<Real> weighted_contrib(n * symm_dofs(3));
  auto f = LAMBDA(LO i) {
    auto x = get_vec<3>(coords, i);
    auto dxc = cross(x - center);
    set_symm(weighted_contrib, i, -masses[i] * (dxc * dxc));
  };
  parallel_for(n, f);
  Real symm[6];
  repro_sum(comm, Reals(weighted_contrib), symm_dofs(3), symm);
  return get_symm<3>(symm, 0);
}

Vector<3> get_axis(CommPtr comm, Reals coords,
    Reals masses, Vector<3> center) {
  auto m = get_matrix(comm, coords, masses, center);
  Matrix<3,3> q;
  Vector<3> l;
  decompose_eigen(m, q, l);
  Int min_i = 0;
  for (Int i = 1; i < 3; ++i) {
    if (l[i] < l[min_i]) {
      min_i = i;
    }
  }
  return q[min_i];
}

Reals get_distances(Reals coords, Vector<3> center, Vector<3> axis) {
  CHECK(coords.size() % 3 == 0);
  auto n = coords.size() / 3;
  Write<Real> distances(n);
  auto f = LAMBDA(LO i) {
    distances[i] = (get_vec<3>(coords, i) - center) * axis;
  };
  parallel_for(n, f);
  return distances;
}

Real get_half_weight(CommPtr comm, Reals masses, Read<I8> marked) {
  auto n = masses.size();
  Write<Real> weighted(n);
  auto f = LAMBDA(LO i) {
    weighted[i] = (Real(marked[i]) * masses[i]);
  };
  parallel_for(n, f);
  return repro_sum(comm, Reals(weighted));
}

Read<I8> mark_half(Reals distances, Real distance) {
  auto n = distances.size();
  Write<I8> marked(n);
  auto f = LAMBDA(LO i) {
    marked[i] = (distances[i] > distance);
  };
  parallel_for(n, f);
  return marked;
}

bool mark_axis_bisection(CommPtr comm, Reals distances,
    Reals masses, Real total_mass, Real tolerance,
    Read<I8>& marked) {
  auto n = distances.size();
  CHECK(n == masses.size());
  auto max_dist = comm->allreduce(max(distances), MAX);
  auto min_dist = comm->allreduce(min(distances), MIN);
  auto range = max2(fabs(min_dist), fabs(max_dist));
  auto step = range / 2.;
  Real distance = 0.;
  for (Int i = 0; i < MANTISSA_BITS; ++i) {
    marked = mark_half(distances, distance);
    auto half_weight = get_half_weight(comm, masses, marked);
    if (are_close(half_weight, total_mass / 2., tolerance, 0.)) {
      return true;
    }
    if (half_weight > total_mass / 2.) {
      distance += step;
    } else {
      distance -= step;
    }
    step /= 2.;
  }
  return false;
}

} //end anonymous namespace

Read<I8> mark_bisection(CommPtr comm,
    Reals coords, Reals masses, Real tolerance) {
  CHECK(coords.size() == masses.size() * 3);
  auto total_mass = repro_sum(comm, masses);
  auto center = get_center(comm, coords, masses, total_mass);
  auto axis = get_axis(comm, coords, masses, center);
  auto dists = get_distances(coords, center, axis);
  Read<I8> marked;
  if (mark_axis_bisection(comm, dists, masses, total_mass, tolerance,
        marked)) {
    return marked;
  }
  // if we couldn't find a decent cutting plane, this may be a highly
  // structured mesh with many points coincident on the cutting plane.
  // perturb the axis vector to escape this scenario
  for (Int i = 0; i < 3 * 2; ++i) {
    auto axis2 = axis;
    axis2[i / 2] += (i % 2) * 1e-6;
    dists = get_distances(coords, center, axis2);
    if (mark_axis_bisection(comm, dists, masses, total_mass, tolerance,
          marked)) {
      return marked;
    }
  }
  std::cerr << "warning: no good inertial bisection found\n";
  return marked;
}

} //end namespace inertia
