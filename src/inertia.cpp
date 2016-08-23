#include "inertia.hpp"

#include <iostream>

#include "access.hpp"
#include "array.hpp"
#include "bipart.hpp"
#include "eigen.hpp"
#include "loop.hpp"
#include "space.hpp"

namespace Omega_h {

namespace inertia {

namespace {

/* since the cross product is involved, we do this all
   in 3D only. */

Vector<3> get_center(
    CommPtr comm, Reals coords, Reals masses, Real total_mass) {
  auto n = masses.size();
  Write<Real> weighted_coords(n * 3);
  auto f = LAMBDA(LO i) {
    set_vector<3>(weighted_coords, i, masses[i] * get_vector<3>(coords, i));
  };
  parallel_for(n, f);
  Vector<3> result;
  repro_sum(comm, Reals(weighted_coords), 3, &result[0]);
  return result / total_mass;
}

Matrix<3, 3> get_matrix(
    CommPtr comm, Reals coords, Reals masses, Vector<3> center) {
  auto n = masses.size();
  Write<Real> weighted_contrib(n * symm_dofs(3));
  auto f = LAMBDA(LO i) {
    auto x = get_vector<3>(coords, i);
    auto dxc = cross(x - center);
    set_symm(weighted_contrib, i, -masses[i] * (dxc * dxc));
  };
  parallel_for(n, f);
  Vector<6> v;
  repro_sum(comm, Reals(weighted_contrib), 6, &v[0]);
  return vector2symm(v);
}

Vector<3> get_axis(CommPtr comm, Reals coords, Reals masses, Vector<3> center) {
  auto m = get_matrix(comm, coords, masses, center);
  auto ed = decompose_eigen(m);
  auto l = ed.l;
  auto q = ed.q;
  Int min_i = 0;
  for (Int i = 1; i < 3; ++i) {
    if (l[i] < l[min_i]) {
      min_i = i;
    }
  }
  return positivize(q[min_i]);
}

Reals get_distances(Reals coords, Vector<3> center, Vector<3> axis) {
  CHECK(coords.size() % 3 == 0);
  auto n = coords.size() / 3;
  Write<Real> distances(n);
  auto f = LAMBDA(LO i) {
    distances[i] = (get_vector<3>(coords, i) - center) * axis;
  };
  parallel_for(n, f);
  return distances;
}

Real get_half_weight(CommPtr comm, Reals masses, Read<I8> marked) {
  auto n = masses.size();
  Write<Real> weighted(n);
  auto f = LAMBDA(LO i) { weighted[i] = (Real(marked[i]) * masses[i]); };
  parallel_for(n, f);
  return repro_sum(comm, Reals(weighted));
}

Read<I8> mark_half(Reals distances, Real distance) {
  auto n = distances.size();
  Write<I8> marked(n);
  auto f = LAMBDA(LO i) { marked[i] = (distances[i] > distance); };
  parallel_for(n, f);
  return marked;
}

bool mark_axis_bisection(CommPtr comm, Reals distances, Reals masses,
    Real total_mass, Real tolerance, Read<I8>& marked) {
  auto n = distances.size();
  CHECK(n == masses.size());
  auto max_dist = comm->allreduce(max(distances), OMEGA_H_MAX);
  auto min_dist = comm->allreduce(min(distances), OMEGA_H_MIN);
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

Read<I8> mark_bisection_internal(CommPtr comm, Reals coords, Reals masses,
    Real tolerance, Vector<3> axis, Vector<3> center, Real total_mass) {
  auto dists = get_distances(coords, center, axis);
  Read<I8> marked;
  if (mark_axis_bisection(comm, dists, masses, total_mass, tolerance, marked)) {
    return marked;
  }
  // if we couldn't find a decent cutting plane, this may be a highly
  // structured mesh with many points coincident on the cutting plane.
  // perturb the axis vector to escape this scenario
  for (Int i = 0; i < 3 * 2; ++i) {
    auto axis2 = axis;
    axis2[i / 2] += (i % 2) ? 1e-3 : -1e-3;
    dists = get_distances(coords, center, axis2);
    if (mark_axis_bisection(
            comm, dists, masses, total_mass, tolerance, marked)) {
      return marked;
    }
  }
  std::cerr << "omega_h warning: no good inertial bisection\n";
  return marked;
}

}  // end anonymous namespace

Read<I8> mark_bisection(
    CommPtr comm, Reals coords, Reals masses, Real tolerance, Vector<3>& axis) {
  CHECK(coords.size() == masses.size() * 3);
  auto total_mass = repro_sum(comm, masses);
  auto center = get_center(comm, coords, masses, total_mass);
  axis = get_axis(comm, coords, masses, center);
  return mark_bisection_internal(
      comm, coords, masses, tolerance, axis, center, total_mass);
}

Read<I8> mark_bisection_given_axis(
    CommPtr comm, Reals coords, Reals masses, Real tolerance, Vector<3> axis) {
  CHECK(coords.size() == masses.size() * 3);
  auto total_mass = repro_sum(comm, masses);
  auto center = get_center(comm, coords, masses, total_mass);
  return mark_bisection_internal(
      comm, coords, masses, tolerance, axis, center, total_mass);
}

Rib recursively_bisect(CommPtr comm, Reals& coords, Reals& masses,
    Remotes& owners, Real tolerance, Rib hints) {
  if (comm->size() == 1) {
    return Rib();
  }
  CHECK(comm->size() % 2 == 0);
  Vector<3> axis;
  Read<I8> marks;
  if (hints.axes.empty()) {
    marks = inertia::mark_bisection(comm, coords, masses, tolerance, axis);
  } else {
    axis = hints.axes.front();
    hints.axes.erase(hints.axes.begin());
    marks = inertia::mark_bisection_given_axis(
        comm, coords, masses, tolerance, axis);
  }
  auto dist = bi_partition(comm, marks);
  coords = dist.exch(coords, 3);
  masses = dist.exch(masses, 1);
  owners = dist.exch(owners, 1);
  auto halfsize = comm->size() / 2;
  comm = comm->split(comm->rank() / halfsize, comm->rank() % halfsize);
  auto out = recursively_bisect(comm, coords, masses, owners, tolerance, hints);
  out.axes.insert(out.axes.begin(), axis);
  return out;
}

}  // end namespace inertia

}  // end namespace Omega_h
