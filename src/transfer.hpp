#ifndef TRANSFER_HPP
#define TRANSFER_HPP

#include "internal.hpp"

namespace Omega_h {

bool is_transfer_ok(XferOpts const& xopts, std::string const& name,
    Omega_h_Xfer type);
bool is_transfer_required(XferOpts const& xopts, std::string const& name,
    Omega_h_Xfer type);

bool should_inherit(Mesh* mesh, XferOpts const& xopts, Int dim, TagBase const* tag);
bool should_interpolate(Mesh* mesh, XferOpts const& xopts, Int dim, TagBase const* tag);
bool should_fit(Mesh* mesh, XferOpts const& opts, Int dim, TagBase const* tag);
bool should_conserve(Mesh* mesh, XferOpts const& opts, Int dim, TagBase const* tag);
bool should_conserve_any(Mesh* mesh, XferOpts const& opts);
bool is_metric(Mesh* mesh, XferOpts const& opts, Int dim, TagBase const* tag);
bool is_size(Mesh* mesh, XferOpts const& opts, Int dim, TagBase const* tag);
bool is_momentum_velocity(Mesh* mesh, XferOpts const& opts, Int dim, TagBase const* tag);
bool has_momentum_velocity(Mesh* mesh, XferOpts const& opts);

void transfer_refine(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh, LOs keys2edges,
    LOs keys2midverts, Int prod_dim, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents);

void transfer_coarsen(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh, LOs keys2verts,
    Adj keys2doms, Int prod_dim, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents);

void transfer_swap(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh, Int prod_dim, LOs keys2edges,
    LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents);

void transfer_copy(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim,
    std::function<bool(TagBase const*)> filter);
void transfer_copy_swap(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh);
void transfer_copy_motion(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh, Int prod_dim);

template <typename T>
void transfer_common3(
    Mesh* new_mesh, Int ent_dim, TagBase const* tagbase, Write<T> new_data);
template <typename T>
void transfer_common2(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs same_ents2old_ents, LOs same_ents2new_ents, TagBase const* tagbase,
    Write<T> new_data);
template <typename T>
void transfer_common(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs same_ents2old_ents, LOs same_ents2new_ents, LOs prods2new_ents,
    TagBase const* tagbase, Read<T> prod_data);
template <typename T>
void transfer_inherit_refine(Mesh* old_mesh, Mesh* new_mesh, LOs keys2edges,
    Int prod_dim, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents, std::string const& name);

void transfer_length(Mesh* old_mesh, Mesh* new_mesh, LOs same_ents2old_ents,
    LOs same_ents2new_ents, LOs prods2new_ents);
void transfer_quality(Mesh* old_mesh, Mesh* new_mesh, LOs same_ents2old_ents,
    LOs same_ents2new_ents, LOs prods2new_ents);
void transfer_pointwise(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents);

#define INST_DECL(T)                                                           \
  extern template void transfer_common3(                                       \
      Mesh* new_mesh, Int ent_dim, TagBase const* tagbase, Write<T> new_data); \
  extern template void transfer_common2(Mesh* old_mesh, Mesh* new_mesh,        \
      Int ent_dim, LOs same_ents2old_ents, LOs same_ents2new_ents,             \
      TagBase const* tagbase, Write<T> new_data);                              \
  extern template void transfer_common(Mesh* old_mesh, Mesh* new_mesh,         \
      Int ent_dim, LOs same_ents2old_ents, LOs same_ents2new_ents,             \
      LOs prods2new_ents, TagBase const* tagbase, Read<T> prod_data);          \
  extern template void transfer_inherit_refine<T>(Mesh * old_mesh,             \
      Mesh * new_mesh, LOs keys2edges, Int prod_dim, LOs keys2prods,           \
      LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents,      \
      std::string const& name);
INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL

}  // end namespace Omega_h

#endif
