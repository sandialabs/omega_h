#include <Omega_h_amr_transfer.hpp>
#include <Omega_h_transfer.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_map.hpp>

namespace Omega_h {

void amr_transfer_linear_interp(Mesh* old_mesh, Mesh* new_mesh,
    Few<LOs, 4> mods2mds, Few<LOs, 4> mods2midverts,
    LOs same_ents2old_ents, LOs same_ents2new_ents,
    TransferOpts opts) {
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tagbase = old_mesh->get_tag(VERT, i);
    if (!should_interpolate(old_mesh, opts, VERT, tagbase)) continue;
    auto ncomps = tagbase->ncomps();
    auto old_data = old_mesh->get_array<Real>(VERT, tagbase->name());
    auto new_data = Write<Real>(new_mesh->nverts() * ncomps);
    for (Int mod_dim = 1; mod_dim <= old_mesh->dim(); ++mod_dim) {
      auto prod_data = average_field(old_mesh, mod_dim, mods2mds[mod_dim], ncomps, old_data);
      map_into(prod_data, mods2midverts[mod_dim], new_data, ncomps);
    }
    transfer_common2(old_mesh, new_mesh, VERT,
        same_ents2old_ents, same_ents2new_ents, tagbase,
        new_data);
  }
}

}
