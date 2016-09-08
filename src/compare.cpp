#include "internal.hpp"

#include <iostream>

#include "algebra.hpp"
#include "array.hpp"
#include "linpart.hpp"
#include "map.hpp"
#include "owners.hpp"
#include "simplices.hpp"

namespace Omega_h {

template <typename T>
struct CompareArrays {
  static bool compare(
      CommPtr comm, Read<T> a, Read<T> b, Real, Real, Int, Int) {
    return comm->reduce_and(a == b);
  }
};

template <>
struct CompareArrays<Real> {
  static bool compare(CommPtr comm, Read<Real> a, Read<Real> b, Real tol,
      Real floor, Int ncomps, Int dim) {
    if (comm->reduce_and(are_close(a, b, tol, floor))) return true;
    /* if floating point arrays are different, we find the value with the
       largest relative difference and print it out for users to determine
       whether this is actually a serious regression
       (and where in the mesh it is most serious)
       or whether tolerances simply need adjusting */
    auto ah = HostRead<Real>(a);
    auto bh = HostRead<Real>(b);
    LO max_i = -1;
    Real max_diff = 0.0;
    for (LO i = 0; i < ah.size(); ++i) {
      auto diff = rel_diff_with_floor(ah[i], bh[i], floor);
      if (diff > max_diff) {
        max_i = i;
        max_diff = diff;
      }
    }
    auto global_start = comm->exscan(GO(ah.size()), OMEGA_H_SUM);
    auto global_max_diff = comm->allreduce(max_diff, OMEGA_H_MAX);
    I32 rank_cand = ArithTraits<I32>::max();
    if (max_diff == global_max_diff) rank_cand = comm->rank();
    auto best_rank = comm->allreduce(rank_cand, OMEGA_H_MIN);
    if (comm->rank() == best_rank) {
      auto global_max_i = global_start + max_i;
      auto ent_global = global_max_i / ncomps;
      auto comp = global_max_i % ncomps;
      std::cout << "max diff at " << singular_names[dim] << " " << ent_global
                << ", comp " << comp << ", values " << ah[max_i] << " vs "
                << bh[max_i] << '\n';
    }
    return false;
  }
};

template <typename T>
static bool compare_copy_data(Int dim, Read<T> a_data, Dist a_dist,
    Read<T> b_data, Dist b_dist, Int ncomps, Real tol, Real floor) {
  auto a_lin_data = reduce_data_to_owners(a_data, a_dist, ncomps);
  auto b_lin_data = reduce_data_to_owners(b_data, b_dist, ncomps);
  CHECK(a_lin_data.size() == b_lin_data.size());
  auto comm = a_dist.parent_comm();
  auto ret = CompareArrays<T>::compare(
      comm, a_lin_data, b_lin_data, tol, floor, ncomps, dim);
  return ret;
}

static Read<GO> get_local_conn(Mesh* mesh, Int dim, bool full) {
  auto low_dim = ((full) ? (dim - 1) : (VERT));
  auto h2l = mesh->ask_down(dim, low_dim);
  auto l_globals = mesh->ask_globals(low_dim);
  auto hl2l_globals = unmap(h2l.ab2b, l_globals, 1);
  return hl2l_globals;
}

Omega_h_Comparison compare_meshes(
    Mesh* a, Mesh* b, Real tol, Real floor, bool verbose, bool full) {
  CHECK(a->comm()->size() == b->comm()->size());
  CHECK(a->comm()->rank() == b->comm()->rank());
  auto comm = a->comm();
  auto should_print = verbose && (comm->rank() == 0);
  if (a->dim() != b->dim()) {
    if (should_print) std::cout << "mesh dimensions differ\n";
    return OMEGA_H_DIFF;
  }
  Omega_h_Comparison result = OMEGA_H_SAME;
  for (Int dim = 0; dim <= a->dim(); ++dim) {
    if (a->nglobal_ents(dim) != b->nglobal_ents(dim)) {
      if (should_print) {
        std::cout << "global " << singular_names[dim] << " counts differ\n";
      }
      return OMEGA_H_DIFF;
    }
    if (!full && (0 < dim) && (dim < a->dim())) continue;
    auto a_globals = a->ask_globals(dim);
    auto b_globals = b->ask_globals(dim);
    auto a_dist = copies_to_linear_owners(comm, a_globals);
    auto b_dist = copies_to_linear_owners(comm, b_globals);
    if (dim > 0) {
      auto a_conn = get_local_conn(a, dim, full);
      auto b_conn = get_local_conn(b, dim, full);
      auto ok = compare_copy_data(
          dim, a_conn, a_dist, b_conn, b_dist, dim + 1, 0.0, 0.0);
      if (!ok) {
        if (should_print) {
          std::cout << singular_names[dim] << " connectivity doesn't match\n";
        }
        result = OMEGA_H_DIFF;
        continue;
      }
    }
    for (Int i = 0; i < a->ntags(dim); ++i) {
      auto tag = a->get_tag(dim, i);
      auto const& name = tag->name();
      if (!b->has_tag(dim, name)) {
        if (should_print) {
          std::cout << singular_names[dim] << " tag \"" << name
                    << "\" exists in first mesh but not second\n";
        }
        result = OMEGA_H_DIFF;
        continue;
      }
      auto ncomps = tag->ncomps();
      bool ok = false;
      switch (tag->type()) {
        case OMEGA_H_I8:
          ok = compare_copy_data(dim, a->get_array<I8>(dim, name), a_dist,
              b->get_array<I8>(dim, name), b_dist, ncomps, tol, floor);
          break;
        case OMEGA_H_I32:
          ok = compare_copy_data(dim, a->get_array<I32>(dim, name), a_dist,
              b->get_array<I32>(dim, name), b_dist, ncomps, tol, floor);
          break;
        case OMEGA_H_I64:
          ok = compare_copy_data(dim, a->get_array<I64>(dim, name), a_dist,
              b->get_array<I64>(dim, name), b_dist, ncomps, tol, floor);
          break;
        case OMEGA_H_F64:
          ok = compare_copy_data(dim, a->get_array<Real>(dim, name), a_dist,
              b->get_array<Real>(dim, name), b_dist, ncomps, tol, floor);
          break;
      }
      if (!ok) {
        if (should_print) {
          std::cout << singular_names[dim] << " tag \"" << name
                    << "\" values are different\n";
        }
        result = OMEGA_H_DIFF;
      }
    }
    for (Int i = 0; i < b->ntags(dim); ++i) {
      auto tag = b->get_tag(dim, i);
      if (!a->has_tag(dim, tag->name())) {
        if (should_print) {
          std::cout << singular_names[dim] << " tag \"" << tag->name()
                    << "\" exists in second mesh but not in first\n";
        }
        if (result == OMEGA_H_SAME) {
          result = OMEGA_H_MORE;
        }
      }
    }
  }
  return result;
}

}  // end namespace Omega_h
