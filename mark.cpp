Read<I8> mark_exposed_sides(Mesh* mesh) {
  auto ns = mesh->nents(mesh->dim() - 1);
  auto s2sc = mesh->ask_up(mesh->dim() - 1, mesh->dim()).a2ab;
  Write<I8> exposed(ns);
  auto f = LAMBDA(LO s) {
    exposed[s] = ((s2sc[s + 1] - s2sc[s]) < 2);
  };
  parallel_for(ns, f);
  return exposed;
}

Read<I8> mark_down(Mesh* mesh, Int high_dim, Int low_dim,
    Read<I8> high_marked) {
  auto l2h = mesh->ask_up(low_dim, high_dim);
  auto l2lh = l2h.a2ab;
  auto lh2h = l2h.ab2b;
  auto nl = mesh->nents(low_dim);
  Write<I8> out(nl, 0);
  auto f = LAMBDA(LO l) {
    for (LO lh = l2lh[l]; lh < l2lh[l + 1]; ++lh)
      if (high_marked[lh2h[lh]])
        out[l] = 1;
  };
  parallel_for(nl, f);
  return mesh->sync_array(low_dim, Read<I8>(out), 1);
}

Read<I8> mark_up(Mesh* mesh, Int low_dim, Int high_dim,
    Read<I8> low_marked) {
  auto l2h = mesh->ask_down(high_dim, low_dim);
  auto deg = simplex_degrees[high_dim][low_dim];
  auto hl2l = l2h.ab2b;
  auto nh = mesh->nents(high_dim);
  Write<I8> out(nh, 0);
  auto f = LAMBDA(LO h) {
    for (Int hhl = 0; hhl < deg; ++hhl) {
      auto l = hl2l[h * deg + hhl];
      if (low_marked[l]) {
        out[h] = 1;
      }
    }
  };
  parallel_for(nh, f);
  return out;
}

Read<I8> mark_by_class_dim(Mesh* mesh, Int ent_dim, Int class_dim) {
  auto e2class_dim = mesh->get_array<I8>(ent_dim, "class_dim");
  return each_eq_to(e2class_dim, static_cast<I8>(class_dim));
}

Read<I8> mark_by_class(Mesh* mesh, Int ent_dim, Int class_dim, I32 class_id) {
  auto e2class_id = mesh->get_array<I32>(ent_dim, "class_id");
  auto id_marks = each_eq_to(e2class_id, class_id);
  return land_each(id_marks, mark_by_class_dim(mesh, ent_dim, class_dim));
}

Read<I8> mark_class_closure(Mesh* mesh, Int ent_dim, Int class_dim, I32 class_id) {
  CHECK(ent_dim <= class_dim);
  auto eq_marks = mark_by_class(mesh, class_dim, class_dim, class_id);
  if (ent_dim == class_dim) return eq_marks;
  return mark_down(mesh, class_dim, ent_dim, eq_marks);
}

Read<I8> mark_class_closures(Mesh* mesh, Int ent_dim,
    std::vector<Int> class_dims, std::vector<I32> class_ids) {
  CHECK(class_dims.size() == class_ids.size());
  auto marks = Read<I8>(mesh->nents(ent_dim), 0);
  for (std::size_t i = 0; i < class_dims.size(); ++i) {
    marks = lor_each(marks, mark_class_closure(
          mesh, ent_dim, class_dims[i], class_ids[i]));
  }
  return marks;
}

Read<I8> mark_dual_layers(Mesh* mesh, Read<I8> marks, Int nlayers) {
  CHECK(mesh->partition() == GHOSTED);
  auto dual = mesh->ask_dual();
  for (Int i = 0; i < nlayers; ++i) {
    marks = graph_reduce(dual, marks, 1, MAX);
    marks = mesh->sync_array(mesh->dim(), marks, 1);
  }
  return marks;
}

GO count_owned_marks(Mesh* mesh, Int ent_dim, Read<I8> marks) {
  if (mesh->could_be_shared(ent_dim)) {
    marks = land_each(marks, mesh->owned(ent_dim));
  }
  return mesh->comm()->allreduce(GO(sum(marks)), SUM);
}

Read<I8> mark_sliver_layers(Mesh* mesh, Real qual_ceil, Int nlayers) {
  CHECK(mesh->partition() == GHOSTED);
  auto quals = mesh->ask_qualities();
  auto elems_are_slivers = each_lt(quals, qual_ceil);
  return mark_dual_layers(mesh, elems_are_slivers, nlayers);
}
