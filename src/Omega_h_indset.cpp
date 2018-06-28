#include "Omega_h_indset.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

namespace indset {

enum {
  NOT_IN = -1,
  UNKNOWN = 0,
  IN = 1,
};

struct Tuples {
  Bytes marks;  // one of IN, NOT_IN, or UNKNOWN
  Reals qualities;
  GOs globals;
};

struct Tuple {
  OMEGA_H_DEVICE Tuple(Tuples const& tuples, LO i) {
    mark = tuples.marks[i];
    quality = tuples.qualities[i];
    global = tuples.globals[i];
  }
  I8 mark;
  Real quality;
  GO global;
  OMEGA_H_INLINE bool operator<(Tuple const& other) {
    if (mark != other.mark) return mark < other.mark;
    if (quality != other.quality) return quality < other.quality;
    return global < other.global;
  }
};

}  // namespace indset

/* Algorithm 5: MIS_parallel
   Nathan Bell, Steven Dalton, and Luke N. Olson.
   "Exposing fine-grained parallelism in algebraic multigrid methods."
   SIAM Journal on Scientific Computing 34.4 (2012): C123-C152.
*/

GOs find_indset(Graph graph, Int distance, Bytes candidates, Reals qualities,
    GOs globals, Dist owners2copies) {
  begin_code("find_indset");
  auto comm = owners2copies.parent_comm();
  auto n = candidates.size();
  auto is_distributed = comm->size() > 1;
  if (is_distributed) OMEGA_H_CHECK(owners2copies.nroots() == n);
  OMEGA_H_CHECK(qualities.size() == n);
  Write<I8> initial_marks(n);
  auto setup = OMEGA_H_LAMBDA(LO i) {
    if (candidates[i])
      initial_marks[i] = indset::UNKNOWN;
    else
      initial_marks[i] = indset::NOT_IN;
  };
  parallel_for(n, setup, "find_indset(setup)");
  auto marks = Bytes(initial_marks);
  Write<GO> owner_globals(n, GO(-1));
  while (get_sum(comm, each_eq_to(marks, I8(indset::UNKNOWN)))) {
    indset::Tuples tuples = {marks, qualities, globals};
    for (Int r = 0; r < distance; ++r) {
      Write<I8> new_marks(n);
      Write<Real> new_qualities(n);
      Write<GO> new_globals(n);
      auto propagate = OMEGA_H_LAMBDA(LO i) {
        auto t = indset::Tuple(tuples, i);
        auto b = graph.a2ab[i];
        auto e = graph.a2ab[i + 1];
        for (auto ij = b; ij < e; ++ij) {
          auto j = graph.ab2b[ij];
          t = max2(t, indset::Tuple(tuples, j));
        }
        new_marks[i] = t.mark;
        new_qualities[i] = t.quality;
        new_globals[i] = t.global;
      };
      parallel_for(n, propagate, "find_indset(propagate)");
      tuples.marks = new_marks;
      tuples.qualities = new_qualities;
      tuples.globals = new_globals;
      if (is_distributed) {
        tuples.marks = owners2copies.exch(tuples.marks, 1);
        tuples.qualities = owners2copies.exch(tuples.qualities, 1);
        tuples.globals = owners2copies.exch(tuples.globals, 1);
      }
    }
    Write<I8> new_marks(n);
    auto accept = OMEGA_H_LAMBDA(LO i) {
      if (marks[i] == indset::UNKNOWN) {
        if (tuples.globals[i] == globals[i]) {
          new_marks[i] = indset::IN;
          owner_globals[i] = globals[i];
        } else if (tuples.marks[i] == indset::IN) {
          new_marks[i] = indset::NOT_IN;
          owner_globals[i] = tuples.globals[i];
        } else {
          new_marks[i] = marks[i];
        }
      } else {
        new_marks[i] = marks[i];
      }
    };
    parallel_for(n, accept, "find_indset(accept)");
    marks = new_marks;
    if (is_distributed) marks = owners2copies.exch(marks, 1);
  }
  end_code();
  return owner_globals;
}

Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Graph graph, Reals qualities, Bytes candidates) {
  auto globals = mesh->globals(ent_dim);
  auto owners2copies = mesh->ask_dist(ent_dim).invert();
  auto distance = 1;
  auto indset_globals = find_indset(
      graph, distance, candidates, qualities, globals, owners2copies);
  return each_eq(indset_globals, globals);
}

Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Reals qualities, Bytes candidates) {
  if (ent_dim == mesh->dim()) return std::move(candidates);
  mesh->owners_have_all_upward(ent_dim);
  OMEGA_H_CHECK(mesh->owners_have_all_upward(ent_dim));
  auto graph = mesh->ask_star(ent_dim);
  return find_indset(mesh, ent_dim, graph, qualities, candidates);
}

}  // end namespace Omega_h
