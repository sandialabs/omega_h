#include <Omega_h_amr_transfer.hpp>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_transfer.hpp>

#include <Omega_h_print.hpp>

namespace Omega_h {

void amr_transfer_linear_interp(Mesh* old_mesh, Mesh* new_mesh,
    Few<LOs, 4> mods2mds, Few<LOs, 4> mods2midverts, LOs same_ents2old_ents,
    LOs same_ents2new_ents, TransferOpts opts) {
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tagbase = old_mesh->get_tag(VERT, i);
    if (!should_interpolate(old_mesh, opts, VERT, tagbase)) continue;
    auto ncomps = tagbase->ncomps();
    auto old_data = old_mesh->get_array<Real>(VERT, tagbase->name());
    auto new_data = Write<Real>(new_mesh->nverts() * ncomps);
    for (Int mod_dim = 1; mod_dim <= old_mesh->dim(); ++mod_dim) {
      auto prod_data =
          average_field(old_mesh, mod_dim, mods2mds[mod_dim], ncomps, old_data);
      map_into(prod_data, mods2midverts[mod_dim], new_data, ncomps);
    }
    transfer_common2(old_mesh, new_mesh, VERT, same_ents2old_ents,
        same_ents2new_ents, tagbase, new_data);
  }
}

static void amr_transfer_levels(Mesh* old_mesh, Mesh* new_mesh,
    Int prod_dim, LOs same_ents2old_ents, LOs same_ents2new_ents,
    Few<LOs, 4> mods2mds, LOs prods2new_ents) {
  auto levels = old_mesh->ask_levels(prod_dim);
  auto tag_data = old_mesh->get_tag<Byte>(prod_dim, "level");
  auto new_data = Write<Byte>(new_mesh->nents(prod_dim), -1);
  Int offset = 0;
  for (Int mod_dim = max2(Int(EDGE), prod_dim); mod_dim <= old_mesh->dim();
      ++mod_dim) {
    auto nprods_per_mod = hypercube_split_degree(mod_dim, prod_dim);
    auto nmods_of_dim = mods2mds[mod_dim].size();
    for (LO i = 0; i < mods2mds[mod_dim].size(); ++i) {
      auto old_idx = mods2mds[mod_dim][i];
      auto old_level = levels[old_idx];
      for (Int j = 0; j < nprods_per_mod; ++j) {
        auto new_idx = prods2new_ents[offset + (i * nprods_per_mod + j)];
        new_data[new_idx] = old_level + 1;
      }
    }
    offset += nprods_per_mod * nmods_of_dim;
  }
  transfer_common2(old_mesh, new_mesh, prod_dim, same_ents2old_ents,
      same_ents2new_ents, tag_data, new_data);
}

void amr_transfer_leaves_and_levels(Mesh* old_mesh, Mesh* new_mesh,
    Int prod_dim, LOs same_ents2old_ents, LOs same_ents2new_ents,
    Few<LOs, 4> mods2mds, LOs prods2new_ents) {
  amr_transfer_levels(old_mesh, new_mesh, prod_dim, same_ents2old_ents,
      same_ents2new_ents, mods2mds, prods2new_ents);
}

}  // namespace Omega_h
