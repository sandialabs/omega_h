#include "internal.hpp"
#include "loop.hpp"
#include "scan.hpp"
#include "simplices.hpp"
#include "transfer_conserve.hpp"

namespace Omega_h {

DEVICE Int add_unique(LO* stack, Int n, LO e, Int stack_max) {
  for (Int i = 0; i < n; ++i)
    if (stack[i] == e) return n;
  CHECK(n < stack_max);
  stack[n++] = e;
  return n;
}

class FindClosureVerts {
  Graph keys2elems;
  LOs elems2verts;
  Int nverts_per_elem;

 public:
  FindClosureVerts(Mesh* mesh, Graph keys2elems_) {
    keys2elems = keys2elems_;
    elems2verts = mesh->ask_elem_verts();
    nverts_per_elem = simplex_degrees[mesh->dim()][VERT];
  }
  LO count() const { return keys2elems.nnodes(); }
  enum { stack_max = (AvgDegree<3, 0, 3>::value + 1) * 2 };
  DEVICE Int run(LO key, LO* stack, Int stack_max) const {
    Int n = 0;
    for (auto key_elem = keys2elems.a2ab[key];
         key_elem < keys2elems.a2ab[key + 1]; ++key_elem) {
      auto elem1 = keys2elems.ab2b[key_elem];
      for (auto elem_vert = elem1 * nverts_per_elem;
           elem_vert < (elem1 + 1) * nverts_per_elem; ++elem_vert) {
        auto vert = elems2verts[elem_vert];
        n = add_unique(stack, n, vert, stack_max);
      }
    }
    return n;
  }
};

template <typename Specialization>
Graph get_graph(Specialization spec) {
  auto n = spec.count();
  auto degrees_w = Write<LO>(n);
  auto count = LAMBDA(LO i) {
    constexpr Int stack_max = Specialization::stack_max;
    LO stack[stack_max];
    degrees_w[i] = spec.run(i, stack, stack_max);
  };
  parallel_for(n, count);
  auto degrees = LOs(degrees_w);
  auto offsets = offset_scan(degrees);
  auto edges = Write<LO>(offsets.last());
  auto fill = LAMBDA(LO i) {
    auto begin = offsets[i];
    auto end = offsets[i + 1];
    auto stack_max = end - begin;
    auto stack = &edges[begin];
    spec.run(i, stack, stack_max);
  };
  parallel_for(n, fill);
  return Graph(offsets, edges);
}

Graph get_closure_verts(Mesh* mesh, Graph keys2elems) {
  FindClosureVerts spec(mesh, keys2elems);
  return get_graph(spec);
}

}  // end namespace Omega_h
