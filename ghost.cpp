Dist get_local_elem_uses2own_verts(Mesh* mesh) {
  auto verts2elems = mesh->ask_up(VERT, mesh->dim());
  auto verts2uses = verts2elems.a2ab;
  auto verts2own_verts = mesh->ask_owners(VERT);
  auto uses2own_verts = expand(verts2own_verts, verts2uses);
  return Dist(mesh->comm(), uses2own_verts, mesh->nverts());
}

Remotes get_local_elem_uses2own_elems(Mesh* mesh) {
  auto verts2elems = mesh->ask_up(VERT, mesh->dim());
  auto uses2elems = verts2elems.ab2b;
  auto elems2own = mesh->ask_owners(mesh->dim());
  return unmap(uses2elems, elems2own);
}

void get_own_verts2own_elem_uses(Mesh* mesh,
    Remotes& serv_uses2own_elems,
    LOs& own_verts2serv_uses) {
  auto local_uses2own_elems = get_local_elem_uses2own_elems(mesh);
  auto local_uses2own_verts = get_local_elem_uses2own_verts(mesh);
  serv_uses2own_elems = local_uses2own_verts.exch(local_uses2own_elems, 1);
  auto own_verts2local_uses = local_uses2own_verts.invert();
  own_verts2serv_uses = own_verts2local_uses.roots2items();
}

Remotes push_elem_uses(
    Remotes serv_uses2own_elems,
    LOs own_verts2serv_uses,
    Dist own_verts2verts) {
  auto nown_verts = own_verts2verts.nroots();
  auto own_verts2serv_verts = own_verts2verts.roots2items();
  auto own_verts2items = multiply_fans(
      own_verts2serv_uses, own_verts2serv_verts);
  auto nitems = own_verts2items.last();
  auto serv_verts2verts = own_verts2verts.items2dests();
  Write<I32> elem_ranks(nitems);
  Write<LO> elem_idxs(nitems);
  Write<I32> vert_ranks(nitems);
  Write<LO> vert_idxs(nitems);
  auto f = LAMBDA(LO ov) {
    auto item = own_verts2items[ov];
    for (auto sv = own_verts2serv_verts[ov];
         sv < own_verts2serv_verts[ov + 1];
         ++sv) {
      for (auto su = own_verts2serv_uses[ov];
           su < own_verts2serv_uses[ov + 1];
           ++su) {
        elem_ranks[item] = serv_uses2own_elems.ranks[su];
        elem_idxs[item] = serv_uses2own_elems.idxs[su];
        vert_ranks[item] = serv_verts2verts.ranks[sv];
        vert_idxs[item] = serv_verts2verts.idxs[sv];
        ++item;
      }
    }
    CHECK(item == own_verts2items[ov + 1]);
  };
  parallel_for(nown_verts, f);
  auto verts2own_verts = own_verts2verts.invert();
  auto nverts = verts2own_verts.nitems();
  auto items2verts_map = Remotes(Read<I32>(vert_ranks), LOs(vert_idxs));
  Dist items2verts(own_verts2verts.parent_comm(),
      items2verts_map, nverts);
  auto items2own_elems = Remotes(Read<I32>(elem_ranks), LOs(elem_idxs));
  return items2verts.exch(items2own_elems, 1);
}

void ghost_mesh(Mesh* mesh) {
  Remotes own_vert_uses2own_elems;
  LOs own_verts2own_vert_uses;
  get_own_verts2own_elem_uses(mesh,
      own_vert_uses2own_elems,
      own_verts2own_vert_uses);
  auto elem_uses = push_elem_uses(
      own_vert_uses2own_elems,
      own_verts2own_vert_uses,
      mesh->ask_dist(VERT).invert());
  auto uses2old_owners = Dist(mesh->comm(), elem_uses, mesh->nelems());
  auto own_elems2elems = find_unique_use_owners(uses2old_owners);
  auto elems2ownners = own_elems2elems.invert();
  Mesh new_mesh;
  migrate_mesh(mesh, new_mesh, elems2ownners, GHOSTED);
  mesh = new_mesh;
}

void partition_by_verts(Mesh* mesh) {
  Remotes own_vert_uses2own_elems;
  LOs own_verts2own_vert_uses;
  get_own_verts2own_elem_uses(mesh,
      own_vert_uses2own_elems,
      own_verts2own_vert_uses);
  auto uses2old_owners = Dist(mesh->comm(),
      own_vert_uses2own_elems, mesh->nelems());
  auto own_elems2elems = find_unique_use_owners(uses2old_owners);
  auto elems2ownners = own_elems2elems.invert();
  Mesh new_mesh;
  migrate_mesh(mesh, new_mesh, elems2ownners, VERTEX_BASED);
  mesh = new_mesh;
}

void partition_by_elems(Mesh* mesh) {
  auto dim = mesh->dim();
  auto all2owners = mesh->ask_owners(dim);
  auto marked_owned = mesh->owned(dim);
  auto owned2all = collect_marked(marked_owned);
  auto owned2owners = unmap(owned2all, all2owners);
  migrate_mesh(mesh, owned2owners);
}
