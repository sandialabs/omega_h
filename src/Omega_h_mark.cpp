#include "Omega_h_mark.hpp"

#include <algorithm>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

Read<I8> mark_exposed_sides(Mesh* mesh) {
  auto ns = mesh->nents(mesh->dim() - 1);
  auto s2sc = mesh->ask_up(mesh->dim() - 1, mesh->dim()).a2ab;
  Write<I8> exposed(ns);
  auto f = OMEGA_H_LAMBDA(LO s) { exposed[s] = ((s2sc[s + 1] - s2sc[s]) < 2); };
  parallel_for(ns, f, "mark_exposed_sides");
  return exposed;
}

Read<I8> mark_down(Graph l2h, Read<I8> high_marked) {
  auto l2lh = l2h.a2ab;
  auto lh2h = l2h.ab2b;
  auto nl = l2lh.size() - 1;
  Write<I8> low_marks_w(nl, 0);
  auto f = OMEGA_H_LAMBDA(LO l) {
    for (LO lh = l2lh[l]; lh < l2lh[l + 1]; ++lh)
      if (high_marked[lh2h[lh]]) low_marks_w[l] = 1;
  };
  parallel_for(nl, std::move(f));
  return low_marks_w;
}

Read<I8> mark_down(
    Mesh* mesh, Int high_dim, Int low_dim, Read<I8> high_marked) {
  OMEGA_H_CHECK(0 <= low_dim);
  OMEGA_H_CHECK(low_dim <= high_dim);
  OMEGA_H_CHECK(high_dim <= 3);
  if (high_dim == low_dim) return high_marked;
  auto l2h = mesh->ask_up(low_dim, high_dim);
  auto low_marks = mark_down(l2h, high_marked);
  if (!mesh->owners_have_all_upward(low_dim)) {
    low_marks = mesh->reduce_array(low_dim, low_marks, 1, OMEGA_H_MAX);
  }
  low_marks = mesh->sync_array(low_dim, low_marks, 1);
  return low_marks;
}

Read<I8> mark_up(Mesh* mesh, Int low_dim, Int high_dim, Read<I8> low_marked) {
  auto l2h = mesh->ask_down(high_dim, low_dim);
  auto deg = element_degree(mesh->family(), high_dim, low_dim);
  auto hl2l = l2h.ab2b;
  auto nh = mesh->nents(high_dim);
  Write<I8> out(nh, 0);
  auto f = OMEGA_H_LAMBDA(LO h) {
    for (Int hhl = 0; hhl < deg; ++hhl) {
      auto l = hl2l[h * deg + hhl];
      if (low_marked[l]) out[h] = 1;
    }
  };
  parallel_for(nh, f, "mark_up");
  return out;
}

Read<I8> mark_adj(Mesh* mesh, Int from_dim, Int to_dim, Read<I8> from_marked) {
  if (from_dim == to_dim) return from_marked;
  if (from_dim < to_dim) return mark_up(mesh, from_dim, to_dim, from_marked);
  if (from_dim > to_dim) return mark_down(mesh, from_dim, to_dim, from_marked);
  OMEGA_H_NORETURN(Read<I8>());
}

Read<I8> mark_up_all(
    Mesh* mesh, Int low_dim, Int high_dim, Read<I8> low_marked) {
  auto l2h = mesh->ask_down(high_dim, low_dim);
  auto deg = element_degree(mesh->family(), high_dim, low_dim);
  auto hl2l = l2h.ab2b;
  auto nh = mesh->nents(high_dim);
  Write<I8> out(nh, 0);
  auto f = OMEGA_H_LAMBDA(LO h) {
    bool all_marked = true;
    for (Int hhl = 0; hhl < deg; ++hhl) {
      auto l = hl2l[h * deg + hhl];
      if (!low_marked[l]) all_marked = false;
    }
    out[h] = I8(all_marked);
  };
  parallel_for(nh, f, "mark_up_all");
  return out;
}

Read<I8> mark_by_class_dim(Mesh* mesh, Int ent_dim, Int class_dim) {
  auto e2class_dim = mesh->get_array<I8>(ent_dim, "class_dim");
  return each_eq_to(e2class_dim, static_cast<I8>(class_dim));
}

Read<I8> mark_by_class(Mesh* mesh, Int ent_dim, Int class_dim, I32 class_id) {
  auto e2class_id = mesh->get_array<ClassId>(ent_dim, "class_id");
  auto id_marks = each_eq_to(e2class_id, class_id);
  return land_each(id_marks, mark_by_class_dim(mesh, ent_dim, class_dim));
}

Read<I8> mark_class_closure(
    Mesh* mesh, Int ent_dim, Int class_dim, ClassId class_id) {
  OMEGA_H_CHECK(ent_dim <= class_dim);
  auto eq_marks = mark_by_class(mesh, class_dim, class_dim, class_id);
  if (ent_dim == class_dim) return eq_marks;
  return mark_down(mesh, class_dim, ent_dim, eq_marks);
}

Read<I8> get_eq_marks(
    Mesh* mesh, Int class_dim, std::vector<ClassId> const& class_ids) {
  auto sorted_class_ids = class_ids;
  std::sort(begin(sorted_class_ids), end(sorted_class_ids));
  HostWrite<LO> h_sorted_class_ids(LO(sorted_class_ids.size()));
  for (size_t i = 0; i < sorted_class_ids.size(); ++i) {
    h_sorted_class_ids[LO(i)] = sorted_class_ids[i];
  }
  auto d_sorted_class_ids = Read<ClassId>(h_sorted_class_ids.write());
  auto nclass_ids = d_sorted_class_ids.size();
  auto eq_class_dims = mesh->get_array<I8>(class_dim, "class_dim");
  auto eq_class_ids = mesh->get_array<LO>(class_dim, "class_id");
  auto neq = mesh->nents(class_dim);
  Write<I8> eq_marks_w(neq);
  auto f = OMEGA_H_LAMBDA(LO eq) {
    eq_marks_w[eq] = I8((eq_class_dims[eq] == I8(class_dim)) &&
                        (-1 != binary_search(d_sorted_class_ids,
                                   eq_class_ids[eq], nclass_ids)));
  };
  parallel_for(neq, f, "mark_class_closures");
  return eq_marks_w;
}

Read<I8> mark_class_closures(Mesh* mesh, Int ent_dim, Int class_dim,
    std::vector<ClassId> const& class_ids) {
  OMEGA_H_CHECK(class_dim >= ent_dim);
  auto eq_marks = get_eq_marks(mesh, class_dim, class_ids);
  auto marks = mark_down(mesh, class_dim, ent_dim, eq_marks);
  return marks;
}

Read<I8> mark_class_closures(Mesh* mesh, Int class_dim,
    std::vector<ClassId> const& class_ids, Graph nodes2ents) {
  OMEGA_H_CHECK(nodes2ents.a2ab.exists());
  OMEGA_H_CHECK(nodes2ents.ab2b.exists());
  auto eq_marks = get_eq_marks(mesh, class_dim, class_ids);
  auto marks = mark_down(nodes2ents, eq_marks);
  return marks;
}

static std::vector<LO> get_dim_class_ids(
    Int class_dim, std::vector<ClassPair> const& class_pairs) {
  std::vector<LO> dim_class_ids;
  for (size_t i = 0; i < class_pairs.size(); ++i) {
    if (class_pairs[i].dim == class_dim) {
      dim_class_ids.push_back(class_pairs[i].id);
    }
  }
  return dim_class_ids;
}

Read<I8> mark_class_closures(
    Mesh* mesh, Int ent_dim, std::vector<ClassPair> const& class_pairs) {
  auto marks = Read<I8>(mesh->nents(ent_dim), I8(0));
  for (Int class_dim = ent_dim; class_dim <= mesh->dim(); ++class_dim) {
    auto dim_class_ids = get_dim_class_ids(class_dim, class_pairs);
    if (dim_class_ids.empty()) continue;
    auto class_dim_marks =
        mark_class_closures(mesh, ent_dim, class_dim, dim_class_ids);
    marks = lor_each(marks, class_dim_marks);
  }
  return marks;
}

Read<I8> mark_class_closures(Mesh* mesh,
    std::vector<ClassPair> const& class_pairs, Graph nodes2ents[4]) {
  auto dim = mesh->dim();
  OMEGA_H_CHECK(nodes2ents[dim].a2ab.exists());
  auto nnodes = nodes2ents[dim].a2ab.size() - 1;
  auto marks = Read<I8>(nnodes, I8(0));
  for (int class_dim = 0; class_dim <= dim; ++class_dim) {
    auto dim_class_ids = get_dim_class_ids(class_dim, class_pairs);
    if (dim_class_ids.empty()) continue;
    auto class_dim_marks = mark_class_closures(
        mesh, class_dim, dim_class_ids, nodes2ents[class_dim]);
    marks = lor_each(marks, class_dim_marks);
  }
  return marks;
}

Read<I8> mark_dual_layers(Mesh* mesh, Read<I8> marks, Int nlayers) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto dual = mesh->ask_dual();
  for (Int i = 0; i < nlayers; ++i) {
    marks = graph_reduce(dual, marks, 1, OMEGA_H_MAX);
    marks = mesh->sync_array(mesh->dim(), marks, 1);
  }
  return marks;
}

GO count_owned_marks(Mesh* mesh, Int ent_dim, Read<I8> marks) {
  if (mesh->could_be_shared(ent_dim)) {
    marks = land_each(marks, mesh->owned(ent_dim));
  }
  return get_sum(mesh->comm(), marks);
}

Read<I8> mark_sliver_layers(Mesh* mesh, Real qual_ceil, Int nlayers) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto quals = mesh->ask_qualities();
  auto elems_are_slivers = each_lt(quals, qual_ceil);
  return mark_dual_layers(mesh, elems_are_slivers, nlayers);
}

}  // end namespace Omega_h
