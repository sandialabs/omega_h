#include "graph.hpp"
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

class FindDistance2Elems {
  Graph keys2elems;
  LOs elems2verts;
  Int nverts_per_elem;
  Graph verts2elems;
  Read<I8> indset;

 public:
  FindDistance2Elems(Mesh* mesh, Int key_dim, Read<I8> unbuffered_indset) {
    keys2elems = mesh->ask_up(key_dim, mesh->dim());
    elems2verts = mesh->ask_elem_verts();
    nverts_per_elem = simplex_degrees[mesh->dim()][VERT];
    verts2elems = mesh->ask_up(VERT, mesh->dim());
    indset = unbuffered_indset;
  }
  LO count() const { return keys2elems.nnodes(); }
  enum { stack_max = 256 };
  DEVICE Int run(LO key, LO* stack, Int stack_max) const {
    if (!indset[key]) return 0;
    Int n = 0;
    for (auto key_elem = keys2elems.a2ab[key];
         key_elem < keys2elems.a2ab[key + 1]; ++key_elem) {
      auto elem1 = keys2elems.ab2b[key_elem];
      for (auto elem_vert = elem1 * nverts_per_elem;
           elem_vert < (elem1 + 1) * nverts_per_elem; ++elem_vert) {
        auto vert = elems2verts[elem_vert];
        for (auto vert_elem = verts2elems.a2ab[vert];
             vert_elem < verts2elems.a2ab[vert + 1]; ++vert_elem) {
          auto elem2 = verts2elems.ab2b[vert_elem];
          n = add_unique(stack, n, elem2, stack_max);
        }
      }
    }
    return n;
  }
};

class FindDistance3Keys {
  Graph keys2buf_elems;
  LOs elems2verts;
  Int nverts_per_elem;
  Graph verts2elems;
  LOs elems2keys;
  Int nkeys_per_elem;
  Read<I8> indset;

 public:
  FindDistance3Keys(Mesh* mesh, Int key_dim, Graph keys2buf_elems_,
      Read<I8> unbuffered_indset) {
    keys2buf_elems = keys2buf_elems_;
    elems2verts = mesh->ask_elem_verts();
    nverts_per_elem = simplex_degrees[mesh->dim()][VERT];
    verts2elems = mesh->ask_up(VERT, mesh->dim());
    elems2keys = mesh->ask_down(mesh->dim(), key_dim).ab2b;
    nkeys_per_elem = simplex_degrees[mesh->dim()][key_dim];
    indset = unbuffered_indset;
  }
  LO count() const { return keys2buf_elems.nnodes(); }
  enum { stack_max = 32 };
  DEVICE Int run(LO key, LO* stack, Int stack_max) const {
    if (!indset[key]) return 0;
    Int n = 0;
    for (auto key_buf_elem = keys2buf_elems.a2ab[key];
         key_buf_elem < keys2buf_elems.a2ab[key + 1]; ++key_buf_elem) {
      auto elem1 = keys2buf_elems.ab2b[key_buf_elem];
      for (auto elem_vert = elem1 * nverts_per_elem;
           elem_vert < (elem1 + 1) * nverts_per_elem; ++elem_vert) {
        auto vert = elems2verts[elem_vert];
        for (auto vert_elem = verts2elems.a2ab[vert];
             vert_elem < verts2elems.a2ab[vert + 1]; ++vert_elem) {
          auto elem2 = verts2elems.ab2b[vert_elem];
          for (auto elem_key = elem2 * nkeys_per_elem;
               elem_key < (elem2 + 1) * nkeys_per_elem; ++elem_key) {
            auto key2 = elems2keys[elem_key];
            if (key2 != key && indset[key2]) {
              n = add_unique(stack, n, key2, stack_max);
            }
          }
        }
      }
    }
    return n;
  }
};

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
  enum { stack_max = (AvgDegree<3,0,3>::value + 1) * 2 };
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

Graph get_buffered_elems(Mesh* mesh, Int key_dim, Read<I8> unbuffered_indset) {
  FindDistance2Elems spec(mesh, key_dim, unbuffered_indset);
  return get_graph(spec);
}

Graph get_buffered_conflicts(
    Mesh* mesh, Int key_dim, Graph keys2buf_elems, Read<I8> unbuffered_indset) {
  FindDistance3Keys spec(mesh, key_dim, keys2buf_elems, unbuffered_indset);
  return get_graph(spec);
}

Graph get_closure_verts(Mesh* mesh, Graph keys2elems) {
  FindClosureVerts spec(mesh, keys2elems);
  return get_graph(spec);
}

Graph get_donor_interior_elems(Mesh* mesh, Int key_dim, LOs keys2kds) {
  auto kds2elems = mesh->ask_up(key_dim, mesh->dim());
  return unmap_graph(keys2kds, kds2elems);
}

Graph get_target_buffer_elems(
    Graph keys2donor_elems, LOs donor_elems2target_elems) {
  auto nedges = keys2donor_elems.nedges();
  Write<I8> keep_w(nedges);
  auto f = LAMBDA(LO edge) {
    keep_w[edge] = (donor_elems2target_elems[keys2donor_elems.ab2b[edge]] >= 0);
  };
  parallel_for(nedges, f);
  auto keep = Read<I8>(keep_w);
  return filter_graph(keys2donor_elems, keep);
}

LOs number_cavity_ents(Mesh* mesh, Graph keys2ents, Int ent_dim) {
  auto nents = mesh->nents(ent_dim);
  auto nkeys = keys2ents.nnodes();
  auto out = Write<LO>(nents, -1);
  auto f = LAMBDA(LO key) {
    Int i = 0;
    for (auto ke = keys2ents.a2ab[key]; ke < keys2ents.a2ab[key + 1]; ++ke) {
      auto ent = keys2ents.ab2b[ke];
      out[ent] = i++;
    }
  };
  parallel_for(nkeys, f);
  return out;
}

}  // end namespace Omega_h
