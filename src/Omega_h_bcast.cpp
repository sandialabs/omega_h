#include "Omega_h_bcast.hpp"

#include "Omega_h_mesh.hpp"

namespace Omega_h {

void bcast_mesh(Mesh* mesh, CommPtr new_comm, bool is_source) {
  if (new_comm->rank() == 0) {
    OMEGA_H_CHECK(is_source);
  }
  I32 family;
  if (is_source) family = mesh->family();
  new_comm->bcast(family);
  if (!is_source) mesh->set_family(Omega_h_Family(family));
  I32 dim;
  if (is_source) dim = mesh->dim();
  new_comm->bcast(dim);
  if (!is_source) mesh->set_dim(dim);
  I32 parting;
  if (is_source) parting = mesh->parting();
  new_comm->bcast(parting);
  if (!is_source) mesh->set_parting(Omega_h_Parting(parting));
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
      if (!is_source) {
        switch (tag_type) {
          case OMEGA_H_I8:
            mesh->add_tag(d, name, ncomps, Read<I8>({}));
            break;
          case OMEGA_H_I32:
            mesh->add_tag(d, name, ncomps, Read<I32>({}));
            break;
          case OMEGA_H_I64:
            mesh->add_tag(d, name, ncomps, Read<I64>({}));
            break;
          case OMEGA_H_F64:
            mesh->add_tag(d, name, ncomps, Read<Real>({}));
            break;
        }
      }
    }
  }
  I32 nsets;
  if (is_source) nsets = I32(mesh->class_sets.size());
  new_comm->bcast(nsets);
  decltype(mesh->class_sets)::iterator it;
  if (is_source) it = mesh->class_sets.begin();
  for (I32 i = 0; i < nsets; ++i) {
    std::string name;
    if (is_source) name = it->first;
    new_comm->bcast_string(name);
    I32 npairs;
    if (is_source) npairs = I32(it->second.size());
    new_comm->bcast(npairs);
    for (I32 j = 0; j < npairs; ++j) {
      ClassPair pair;
      if (is_source) pair = it->second[std::size_t(j)];
      new_comm->bcast(pair.dim);
      new_comm->bcast(pair.id);
      if (!is_source) mesh->class_sets[name].push_back(pair);
    }
    if (is_source) ++it;
  }
}

}  // end namespace Omega_h
