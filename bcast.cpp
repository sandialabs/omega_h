void bcast_mesh(Mesh* mesh, CommPtr new_comm, bool is_source) {
  if (new_comm->rank() == 0) {
    CHECK(is_source);
  }
  I32 dim;
  if (is_source) dim = mesh->dim();
  new_comm->bcast(dim);
  if (!is_source) mesh->set_dim(dim);
  I32 partition;
  if (is_source) partition = mesh->partition();
  new_comm->bcast(partition);
  if (!is_source) mesh->set_partition(static_cast<Partition>(partition));
  I32 keep_canon;
  if (is_source) keep_canon = mesh->keeps_canonical_globals();
  new_comm->bcast(keep_canon);
  if (!is_source) mesh->keep_canonical_globals(keep_canon);
  if (!is_source) mesh->set_verts(0);
  for (Int d = 0; d <= dim; ++d) {
    if (d > VERT && !is_source) {
      Adj down(LOs({}));
      if (d - 1 != VERT) down.codes = Read<I8>({});
      mesh->set_ents(d, down);
    }
    I32 ntags;
    if (is_source) ntags = mesh->ntags(d);
    new_comm->bcast(ntags);
    for (Int i = 0; i < ntags; ++i) {
      TagBase const* tag = nullptr;
      if (is_source) tag = mesh->get_tag(d, i);
      std::string name;
      if (is_source) name = tag->name();
      new_comm->bcast_string(name);
      I32 ncomps;
      if (is_source) ncomps = tag->ncomps();
      new_comm->bcast(ncomps);
      I32 tag_type;
      if (is_source) tag_type = tag->type();
      new_comm->bcast(tag_type);
      I32 tag_xfer;
      if (is_source) tag_xfer = tag->xfer();
      new_comm->bcast(tag_xfer);
      osh_xfer xfer = static_cast<osh_xfer>(tag_xfer);
      if (!is_source) {
        switch (tag_type) {
          case OSH_I8: mesh->add_tag(d, name, ncomps, xfer, Read<I8>({}));
            break;
          case OSH_I32: mesh->add_tag(d, name, ncomps, xfer, Read<I32>({}));
            break;
          case OSH_I64: mesh->add_tag(d, name, ncomps, xfer, Read<I64>({}));
            break;
          case OSH_F64: mesh->add_tag(d, name, ncomps, xfer, Read<Real>({}));
            break;
        }
      }
    }
  }
}
