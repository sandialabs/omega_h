template <typename T>
static void transfer_common(
    Mesh& old_mesh,
    Mesh& new_mesh,
    Int ent_dim,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs prods2new_ents,
    TagBase const* tagbase,
    Read<T> prod_data) {
  auto name = tagbase->name();
  auto ncomps = tagbase->ncomps();
  auto xfer = tagbase->xfer();
  auto old_data = old_mesh.get_array<T>(ent_dim, name);
  auto same_data = unmap(same_ents2old_ents, old_data, ncomps);
  auto nnew_ents = new_mesh.nents(ent_dim);
  auto new_data = Write<T>(nnew_ents * ncomps);
  map_into(same_data, same_ents2new_ents, new_data, ncomps);
  map_into(prod_data, prods2new_ents, new_data, ncomps);
  new_mesh.add_tag(ent_dim, name, ncomps, xfer, Read<T>(new_data));
}

void transfer_linear_interp(Mesh& old_mesh, Mesh& new_mesh,
    LOs keys2edges,
    LOs keys2midverts,
    LOs same_verts2old_verts,
    LOs same_verts2new_verts) {
  for (Int i = 0; i < old_mesh.ntags(VERT); ++i) {
    auto tagbase = old_mesh.get_tag(VERT, i);
    if (tagbase->xfer() == OSH_LINEAR_INTERP) {
      auto ncomps = tagbase->ncomps();
      auto old_data = old_mesh.get_array<Real>(VERT, tagbase->name());
      auto prod_data = average_field(old_mesh, EDGE, keys2edges, ncomps, old_data);
      transfer_common(old_mesh, new_mesh, VERT,
          same_verts2old_verts, same_verts2new_verts,
          keys2midverts, tagbase, prod_data);
    }
  }
}

void transfer_metric(Mesh& old_mesh, Mesh& new_mesh,
    LOs keys2edges,
    LOs keys2midverts,
    LOs same_verts2old_verts,
    LOs same_verts2new_verts) {
  for (Int i = 0; i < old_mesh.ntags(VERT); ++i) {
    auto tagbase = old_mesh.get_tag(VERT, i);
    if (tagbase->xfer() == OSH_METRIC) {
      auto old_data = old_mesh.get_array<Real>(VERT, tagbase->name());
      auto prod_data = average_metric(old_mesh, EDGE, keys2edges, old_data);
      transfer_common(old_mesh, new_mesh, VERT,
          same_verts2old_verts, same_verts2new_verts,
          keys2midverts, tagbase, prod_data);
    }
  }
}

template <typename T>
static void transfer_inherit_refine(Mesh& old_mesh, Mesh& new_mesh,
    LOs keys2edges,
    Int prod_dim,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    std::string const& name) {
  auto old_tag = old_mesh.get_tag<T>(prod_dim, name);
  auto ncomps = old_tag->ncomps();
  auto nprods = keys2prods.last();
  auto prod_data = Write<T>(nprods * ncomps);
  auto nkeys = keys2edges.size();
  /* transfer pairs */
  if (prod_dim > VERT) {
    auto dom_dim = prod_dim;
    auto dom_data = old_mesh.get_array<T>(dom_dim, name);
    auto edges2doms = old_mesh.ask_graph(EDGE, dom_dim);
    auto edges2edge_doms = edges2doms.a2ab;
    auto edge_doms2doms = edges2doms.ab2b;
    auto f = LAMBDA(LO key) {
      auto edge = keys2edges[key];
      auto prod = keys2prods[key];
      for (auto edge_dom = edges2edge_doms[edge];
           edge_dom < edges2edge_doms[edge + 1];
           ++edge_dom) {
        auto dom = edge_doms2doms[edge_dom];
        for (Int pair = 0; pair < 2; ++pair) {
          for (Int comp = 0; comp < ncomps; ++comp) {
            prod_data[prod * ncomps + comp] =
              dom_data[dom * ncomps + comp];
          }
          ++prod;
        }
      }
    };
    parallel_for(nkeys, f);
  }
  if (prod_dim < old_mesh.dim()) {
    auto dom_dim = prod_dim + 1;
    auto dom_data = old_mesh.get_array<T>(dom_dim, name);
    auto edges2doms = old_mesh.ask_graph(EDGE, dom_dim);
    auto edges2edge_doms = edges2doms.a2ab;
    auto edge_doms2doms = edges2doms.ab2b;
    auto f = LAMBDA(LO key) {
      auto edge = keys2edges[key];
      auto ndoms = edges2edge_doms[edge + 1] - edges2edge_doms[edge];
      auto prod = keys2prods[key + 1] - ndoms;
      for (auto edge_dom = edges2edge_doms[edge];
           edge_dom < edges2edge_doms[edge + 1];
           ++edge_dom) {
        auto dom = edge_doms2doms[edge_dom];
        for (Int comp = 0; comp < ncomps; ++comp) {
          prod_data[prod * ncomps + comp] =
            dom_data[dom * ncomps + comp];
        }
        ++prod;
      }
    };
    parallel_for(nkeys, f);
  }
  transfer_common(old_mesh, new_mesh, prod_dim,
      same_ents2old_ents, same_ents2new_ents, prods2new_ents,
      old_tag, Read<T>(prod_data));
}

void transfer_inherit_refine(Mesh& old_mesh, Mesh& new_mesh,
    LOs keys2edges,
    Int prod_dim,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  for (Int i = 0; i < old_mesh.ntags(prod_dim); ++i) {
    auto tagbase = old_mesh.get_tag(prod_dim, i);
    if (tagbase->xfer() == OSH_INHERIT) {
      switch(tagbase->type()) {
      case OSH_I8:
        transfer_inherit_refine<I8>(old_mesh, new_mesh,
            keys2edges, prod_dim, keys2prods, prods2new_ents,
            same_ents2old_ents, same_ents2new_ents, tagbase->name());
        break;
      case OSH_I32:
        transfer_inherit_refine<I32>(old_mesh, new_mesh,
            keys2edges, prod_dim, keys2prods, prods2new_ents,
            same_ents2old_ents, same_ents2new_ents, tagbase->name());
        break;
      case OSH_I64:
        transfer_inherit_refine<I64>(old_mesh, new_mesh,
            keys2edges, prod_dim, keys2prods, prods2new_ents,
            same_ents2old_ents, same_ents2new_ents, tagbase->name());
        break;
      case OSH_F64:
        transfer_inherit_refine<Real>(old_mesh, new_mesh,
            keys2edges, prod_dim, keys2prods, prods2new_ents,
            same_ents2old_ents, same_ents2new_ents, tagbase->name());
        break;
      }
    }
  }
}

void transfer_length(Mesh& old_mesh, Mesh& new_mesh,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs prods2new_ents) {
  for (Int i = 0; i < old_mesh.ntags(EDGE); ++i) {
    auto tagbase = old_mesh.get_tag(EDGE, i);
    if (tagbase->xfer() == OSH_LENGTH) {
      auto prod_data = measure_edges(new_mesh, prods2new_ents);
      transfer_common(old_mesh, new_mesh, EDGE,
          same_ents2old_ents, same_ents2new_ents, prods2new_ents,
          tagbase, prod_data);
    }
  }
}

void transfer_quality(Mesh& old_mesh, Mesh& new_mesh,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs prods2new_ents) {
  auto dim = old_mesh.dim();
  for (Int i = 0; i < old_mesh.ntags(dim); ++i) {
    auto tagbase = old_mesh.get_tag(dim, i);
    if (tagbase->xfer() == OSH_QUALITY) {
      auto prod_data = measure_qualities(new_mesh, prods2new_ents);
      transfer_common(old_mesh, new_mesh, dim,
          same_ents2old_ents, same_ents2new_ents, prods2new_ents,
          tagbase, prod_data);
    }
  }
}

void transfer_refine(Mesh& old_mesh, Mesh& new_mesh,
    LOs keys2edges,
    LOs keys2midverts,
    Int prod_dim,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  transfer_inherit_refine(old_mesh, new_mesh, keys2edges, prod_dim,
      keys2prods, prods2new_ents,
      same_ents2old_ents, same_ents2new_ents);
  if (prod_dim == VERT) {
    transfer_linear_interp(old_mesh, new_mesh, keys2edges, keys2midverts,
        same_ents2old_ents, same_ents2new_ents);
    transfer_metric(old_mesh, new_mesh, keys2edges, keys2midverts,
        same_ents2old_ents, same_ents2new_ents);
  } else if (prod_dim == EDGE) {
    transfer_length(old_mesh, new_mesh,
        same_ents2old_ents, same_ents2new_ents, prods2new_ents);
  } else if (prod_dim == old_mesh.dim()) {
    transfer_quality(old_mesh, new_mesh,
        same_ents2old_ents, same_ents2new_ents, prods2new_ents);
  }
}
