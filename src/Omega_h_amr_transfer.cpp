#include <Omega_h_amr.hpp>
#include <Omega_h_amr_transfer.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_transfer.hpp>

namespace Omega_h {

namespace amr {

void transfer_linear_interp(Mesh* old_mesh, Mesh* new_mesh,
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

void transfer_levels(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim,
    Few<LOs, 4> mods2mds, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  old_mesh->ask_levels(prod_dim);
  auto dim = old_mesh->dim();
  auto new_data = Write<Byte>(new_mesh->nents(prod_dim), -1);
  Int offset = 0;
  for (Int mod_dim = max2(Int(EDGE), prod_dim); mod_dim <= dim; ++mod_dim) {
    auto old_mod_data = old_mesh->ask_levels(mod_dim);
    auto nprods_per_mod = hypercube_split_degree(mod_dim, prod_dim);
    auto nmods_of_dim = mods2mds[mod_dim].size();
    auto f = OMEGA_H_LAMBDA(LO mod) {
      auto old_idx = mods2mds[mod_dim][mod];
      auto old_level = old_mod_data[old_idx];
      for (Int prod = 0; prod < nprods_per_mod; ++prod) {
        auto new_idx = prods2new_ents[offset + (mod * nprods_per_mod + prod)];
        new_data[new_idx] = old_level + 1;
      }
    };
    parallel_for(nmods_of_dim, f);
    offset += nprods_per_mod * nmods_of_dim;
  }
  auto tag = old_mesh->get_tag<Byte>(prod_dim, "level");
  transfer_common2(old_mesh, new_mesh, prod_dim, same_ents2old_ents,
      same_ents2new_ents, tag, new_data);
}

void transfer_leaves(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim,
    Few<LOs, 4> mods2mds, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents, LOs old_ents2new_ents) {
  if (prod_dim == VERT) return;
  auto dim = old_mesh->dim();
  auto ncomps = 1;
  auto old_data = old_mesh->ask_leaves(prod_dim);
  auto new_data = Write<Byte>(new_mesh->nents(prod_dim), -1);
  auto same_data = read(unmap(same_ents2old_ents, old_data, ncomps));
  map_into(same_data, same_ents2new_ents, new_data, ncomps);
  Int offset = 0;
  for (Int mod_dim = max2(Int(EDGE), prod_dim); mod_dim <= dim; ++mod_dim) {
    auto nprods_per_mod = hypercube_split_degree(mod_dim, prod_dim);
    auto nmods_of_dim = mods2mds[mod_dim].size();
    auto f = OMEGA_H_LAMBDA(LO mod) {
      for (Int prod = 0; prod < nprods_per_mod; ++prod) {
        auto new_idx = prods2new_ents[offset + (mod * nprods_per_mod + prod)];
        new_data[new_idx] = 1;
      }
    };
    parallel_for(nmods_of_dim, f);
    offset += nprods_per_mod * nmods_of_dim;
  }
  auto mods2new_ents = unmap(mods2mds[prod_dim], old_ents2new_ents, ncomps);
  auto f = OMEGA_H_LAMBDA(LO mod) {
    auto mod_idx = mods2new_ents[mod];
    new_data[mod_idx] = 0;
  };
  parallel_for(mods2mds[prod_dim].size(), f);
  new_mesh->add_tag<Byte>(prod_dim, "leaf", 1, new_data, true);
}

void transfer_parents(Mesh* old_mesh, Mesh* new_mesh, Few<LOs, 4> mods2mds,
    Few<LOs, 4> prods2new_ents, Few<LOs, 4> same_ents2old_ents,
    Few<LOs, 4> same_ents2new_ents, Few<LOs, 4> old_ents2new_ents) {
  auto dim = old_mesh->dim();
  for (Int prod_dim = 0; prod_dim <= dim; ++prod_dim) {
    auto old_p_data = (old_mesh->ask_parents(prod_dim)).parent_idx;
    auto old_c_data = (old_mesh->ask_parents(prod_dim)).codes;
    auto new_p_data = Write<LO>(new_mesh->nents(prod_dim), -1);
    auto new_c_data = Write<I8>(new_mesh->nents(prod_dim), 0);
    auto same_p_data = unmap(same_ents2old_ents[prod_dim], old_p_data, 1);
    auto same_c_data = read(unmap(same_ents2old_ents[prod_dim], old_c_data, 1));
    auto lambda1 = OMEGA_H_LAMBDA(LO same_ent) {
      auto old_parent = same_p_data[same_ent];
      if (old_parent < 0) return;
      auto code = same_c_data[same_ent];
      auto parent_dim = code_parent_dim(code);
      auto new_parent = old_ents2new_ents[parent_dim][old_parent];
      same_p_data[same_ent] = new_parent;
    };
    parallel_for(same_p_data.size(), lambda1);
    map_into(read(same_p_data), same_ents2new_ents[prod_dim], new_p_data, 1);
    map_into(same_c_data, same_ents2new_ents[prod_dim], new_c_data, 1);
    Int offset = 0;
    for (Int mod_dim = max2(Int(EDGE), prod_dim); mod_dim <= dim; ++mod_dim) {
      auto mods2new_ents =
          unmap(mods2mds[mod_dim], old_ents2new_ents[mod_dim], 1);
      auto nprods_per_mod = hypercube_split_degree(mod_dim, prod_dim);
      auto nmods_of_dim = mods2mds[mod_dim].size();
      auto lambda2 = OMEGA_H_LAMBDA(LO mod) {
        for (Int prod = 0; prod < nprods_per_mod; ++prod) {
          auto idx = offset + (mod * nprods_per_mod + prod);
          auto prod_idx = prods2new_ents[prod_dim][idx];
          new_p_data[prod_idx] = mods2new_ents[mod];
          new_c_data[prod_idx] = amr::make_code(prod, mod_dim);
        }
      };
      parallel_for(nmods_of_dim, lambda2);
      offset += nprods_per_mod * nmods_of_dim;
    }
    Parents parents(new_p_data, new_c_data);
    new_mesh->set_parents(prod_dim, parents);
  }
}

template <typename T>
void transfer_inherit(Mesh* old_mesh, Mesh* new_mesh,
    Few<LOs, 4> const prods2new_ents, Few<LOs, 4> same_ents2old_ents,
    Few<LOs, 4> same_ents2new_ents, std::string const& name) {
  Few<Write<T>, 4> new_data;
  auto ncomps = old_mesh->get_tag<T>(0, name)->ncomps();
  for (Int prod_dim = 0; prod_dim <= old_mesh->dim(); ++prod_dim) {
    new_data[prod_dim] = Write<T>(new_mesh->nents(prod_dim), -1);
    auto old_data = old_mesh->get_array<T>(prod_dim, name);
    auto same_data =
        read(unmap(same_ents2old_ents[prod_dim], old_data, ncomps));
    map_into(
        same_data, same_ents2new_ents[prod_dim], new_data[prod_dim], ncomps);
  }
  for (Int prod_dim = 0; prod_dim <= old_mesh->dim(); ++prod_dim) {
    auto parents = new_mesh->ask_parents(prod_dim);
    auto f = OMEGA_H_LAMBDA(LO prod) {
      auto prod_ent = prods2new_ents[prod_dim][prod];
      auto parent_ent = parents.parent_idx[prod_ent];
      auto code = parents.codes[prod_ent];
      auto parent_dim = code_parent_dim(code);
      new_data[prod_dim][prod_ent] = new_data[parent_dim][parent_ent];
    };
    parallel_for(prods2new_ents[prod_dim].size(), f);
    new_mesh->add_tag(
        prod_dim, name, ncomps, Read<T>(new_data[prod_dim]), true);
  }
}

static void transfer_inherit(Mesh* old_mesh, Mesh* new_mesh,
    Few<LOs, 4> prods2new_ents, Few<LOs, 4> same_ents2old_ents,
    Few<LOs, 4> same_ents2new_ents, TagBase const* tagbase) {
  auto name = tagbase->name();
  switch (tagbase->type()) {
    case OMEGA_H_I8:
      amr::transfer_inherit<I8>(old_mesh, new_mesh, prods2new_ents,
          same_ents2old_ents, same_ents2new_ents, name);
      break;
    case OMEGA_H_I32:
      amr::transfer_inherit<I32>(old_mesh, new_mesh, prods2new_ents,
          same_ents2old_ents, same_ents2new_ents, name);
      break;
    case OMEGA_H_I64:
      amr::transfer_inherit<I64>(old_mesh, new_mesh, prods2new_ents,
          same_ents2old_ents, same_ents2new_ents, name);
      break;
    case OMEGA_H_F64:
      amr::transfer_inherit<Real>(old_mesh, new_mesh, prods2new_ents,
          same_ents2old_ents, same_ents2new_ents, name);
      break;
  }
}

/* inherited quantities come directly from parents, which have
   multiple dimensions, so we enforce that inherited tags exist
   across all valid dimensions */
static void validate_tag(Mesh* m, TagBase const* tagbase) {
  auto name = tagbase->name();
  for (Int d = 1; d <= m->dim(); ++d) {
    OMEGA_H_CHECK(m->has_tag(d, name));
    auto dtagbase = m->get_tagbase(d, name);
    OMEGA_H_CHECK(tagbase->type() == dtagbase->type());
    OMEGA_H_CHECK(tagbase->ncomps() == dtagbase->ncomps());
  }
}

void transfer_inherit(Mesh* old_mesh, Mesh* new_mesh,
    Few<LOs, 4> prods2new_ents, Few<LOs, 4> same_ents2old_ents,
    Few<LOs, 4> same_ents2new_ents, TransferOpts const& opts) {
  for (Int i = 0; i < old_mesh->ntags(0); ++i) {
    auto tagbase = old_mesh->get_tag(0, i);
    if (should_inherit(old_mesh, opts, 0, tagbase)) {
      validate_tag(old_mesh, tagbase);
      amr::transfer_inherit(old_mesh, new_mesh, prods2new_ents,
          same_ents2old_ents, same_ents2new_ents, tagbase);
    }
  }
}

}  // namespace amr

}  // namespace Omega_h
