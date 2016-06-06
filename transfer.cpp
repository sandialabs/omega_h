void transfer_linear_interp(Mesh& old_mesh, Mesh& new_mesh,
    LOs key_verts2verts,
    LOs keys2midverts,
    LOs same_verts2new_verts) {
  for (Int i = 0; i < old_mesh.ntags(VERT); ++i) {
    auto tagbase = old_mesh.get_tag(VERT, i);
    if (tagbase->type() == OSH_F64 && tagbase->xfer() == OSH_LINEAR_INTERP) {
      auto name = tagbase->name();
      auto ncomps = tagbase->ncomps();
      auto old_data = old_mesh.get_array<Real>(VERT, name);
      auto same_data = old_data;
      auto nnew_verts = new_mesh.nents(VERT);
      auto new_data = Write<Real>(nnew_verts * ncomps);
      map_into(same_data, same_verts2new_verts, new_data, ncomps);
      auto prod_data = average_field(2, key_verts2verts, ncomps, old_data);
      map_into(prod_data, keys2midverts, new_data, ncomps);
      new_mesh.add_tag(VERT, name, ncomps, OSH_LINEAR_INTERP, Reals(new_data));
    }
  }
}
