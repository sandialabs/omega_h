template <typename T>
static void transfer_common(
    Mesh& old_mesh,
    Mesh& new_mesh,
    Int ent_dim,
    LOs same_ents2new_ents,
    LOs prods2new_ents,
    TagBase const* tagbase,
    Read<T> prod_data) {
  auto name = tagbase->name();
  auto ncomps = tagbase->ncomps();
  auto xfer = tagbase->xfer();
  auto old_data = old_mesh.get_array<T>(ent_dim, name);
  auto same_data = old_data;
  auto nnew_ents = new_mesh.nents(ent_dim);
  auto new_data = Write<T>(nnew_ents * ncomps);
  map_into(same_data, same_ents2new_ents, new_data, ncomps);
  map_into(prod_data, prods2new_ents, new_data, ncomps);
  new_mesh.add_tag(ent_dim, name, ncomps, xfer, Read<T>(new_data));
}

void transfer_linear_interp(Mesh& old_mesh, Mesh& new_mesh,
    LOs keys2edges,
    LOs keys2midverts,
    LOs same_verts2new_verts) {
  for (Int i = 0; i < old_mesh.ntags(VERT); ++i) {
    auto tagbase = old_mesh.get_tag(VERT, i);
    if (tagbase->xfer() == OSH_LINEAR_INTERP) {
      auto ncomps = tagbase->ncomps();
      auto old_data = old_mesh.get_array<Real>(VERT, tagbase->name());
      auto prod_data = average_field(old_mesh, EDGE, keys2edges, ncomps, old_data);
      transfer_common(old_mesh, new_mesh, VERT,
          same_verts2new_verts, keys2midverts, tagbase, prod_data);
    }
  }
}

void transfer_metric(Mesh& old_mesh, Mesh& new_mesh,
    LOs keys2edges,
    LOs keys2midverts,
    LOs same_verts2new_verts) {
  for (Int i = 0; i < old_mesh.ntags(VERT); ++i) {
    auto tagbase = old_mesh.get_tag(VERT, i);
    if (tagbase->xfer() == OSH_METRIC) {
      auto old_data = old_mesh.get_array<Real>(VERT, tagbase->name());
      auto prod_data = average_metric(old_mesh, EDGE, keys2edges, old_data);
      transfer_common(old_mesh, new_mesh, VERT,
          same_verts2new_verts, keys2midverts, tagbase, prod_data);
    }
  }
}
