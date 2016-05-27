void bcast_mesh(Mesh& mesh, CommPtr new_comm) {
  CommPtr old_comm = mesh.comm();
  bool is_root = ((old_comm) && (old_comm->rank() == 0));
  if (is_root)
    CHECK(new_comm->rank() == 0);
  I32 dim;
  if (is_root)
    dim = mesh.dim();
  new_comm->bcast(dim);
  if (!is_root)
    mesh.set_dim(dim);
  for (Int d = 0; d <= dim; ++d) {
    I32 ntags;
    if (is_root)
      ntags = mesh.ntags(d);
    new_comm->bcast(ntags);
    for (Int i = 0; i < ntags; ++i) {
      TagBase const* tag = nullptr;
      if (is_root)
        tag = mesh.get_tag(d, i);
      std::string name;
      if (is_root)
        name = tag->name();
      new_comm->bcast_string(name);
      I32 ncomps;
      if (is_root)
        ncomps = tag->ncomps();
      new_comm->bcast(ncomps);
      I32 tag_type;
      if (is_root)
        tag_type = tag->type();
      new_comm->bcast(tag_type);
      if (!is_root) {
        switch (tag_type) {
          case OSH_I8: mesh.add_tag<I8>(d, name, ncomps); break;
          case OSH_I32: mesh.add_tag<I32>(d, name, ncomps); break;
          case OSH_I64: mesh.add_tag<I64>(d, name, ncomps); break;
          case OSH_F64: mesh.add_tag<Real>(d, name, ncomps); break;
        }
      }
    }
  }
}
