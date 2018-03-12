#include <Omega_h_amr.hpp>
#include <Omega_h_amr_topology.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_modify.hpp>

namespace Omega_h {

static void amr_refine_ghosted(Mesh* mesh) {
  Few<LOs, 4> mods2mds;
  for (Int mod_dim = 0; mod_dim <= mesh->dim(); ++mod_dim) {
    auto mds_are_mods = mesh->get_array<Byte>(mod_dim, "refine");
    mods2mds[mod_dim] = collect_marked(mds_are_mods);
  }
  for (Int prod_dim = 0; prod_dim <= mesh->dim(); ++prod_dim) {
    Few<LOs, 4> mods2nprods;
    Few<bool, 4> mods_have_prods;
    mods_have_prods[0] = false;
    for (Int mod_dim = 1; mod_dim <= mesh->dim(); ++mod_dim) {
      auto nprods_per_mods = hypercube_split_degree(mod_dim, prod_dim);
      mods2nprods[mod_dim] = LOs(mods2mds[mod_dim].size(), nprods_per_mods);
      mods_have_prods[mod_dim] = true;
    }
    auto rep2md_orders = get_rep2md_order(mesh, prod_dim, mods2mds, mods2nprods, mods_have_prods);
    auto name = std::string("rep_") + hypercube_singular_name(prod_dim) + "2md_order";
    for (Int mod_dim = 1; mod_dim <= mesh->dim(); ++mod_dim) {
      mesh->add_tag(mod_dim, name, 1, rep2md_orders[mod_dim]);
    }
  }
}

static void amr_refine_elem_based(Mesh* mesh) {
  auto prod_counts = count_amr(mesh);
  Few<Bytes, 4> mds_are_mods;
  Few<LOs, 4> mods2mds;
  Few<LOs, 4> mds2mods;
  for (Int mod_dim = 0; mod_dim <= mesh->dim(); ++mod_dim) {
    mds_are_mods[mod_dim] = mesh->get_array<Byte>(mod_dim, "refine");
    mods2mds[mod_dim] = collect_marked(mds_are_mods[mod_dim]);
    mds2mods[mod_dim] = invert_injective_map(mods2mds[mod_dim], mesh->nents(mod_dim));
  }
  auto new_mesh = mesh->copy_meta();
  Few<LOs, 4> mods2midverts;
  LOs old_lows2new_lows;
  for (Int prod_dim = 0; prod_dim <= mesh->dim(); ++prod_dim) {
    LOs prods2verts;
    if (prod_dim != VERT) {
      prods2verts = get_amr_topology(mesh, prod_dim, prod_counts[prod_dim], mods2mds, mds2mods, mods2midverts);
    }
    Few<LOs, 4> mods2prods;
    for (Int mod_dim = max2(Int(EDGE), prod_dim); mod_dim <= mesh->dim(); ++mod_dim) {
      auto nprods_per_mod = hypercube_split_degree(mod_dim, prod_dim);
      mods2prods[mod_dim] = LOs(mods2mds[mod_dim].size(), nprods_per_mod);
    }
    LOs prods2new_ents;
    LOs same_ents2old_ents;
    LOs same_ents2new_ents;
    LOs old_ents2new_ents;
    modify_ents(mesh, &new_mesh, prod_dim, mods2mds, mds_are_mods, mods2prods, prods2verts, old_lows2new_lows,
        /*keep_mods*/true, /*mods_can_be_shared*/true, &prods2new_ents, &same_ents2old_ents,
        &same_ents2new_ents, &old_ents2new_ents);
    /*TODO: transfer!*/
    if (prod_dim == VERT) {
      LO offset = 0;
      for (Int mod_dim = EDGE; mod_dim <= mesh->dim(); ++mod_dim) {
        OMEGA_H_CHECK(hypercube_split_degree(mod_dim, prod_dim) == 1);
        auto nmods_of_dim = mods2mds[mod_dim].size();
        auto begin = offset;
        auto end = begin + nmods_of_dim;
        mods2midverts[mod_dim] = unmap_range(begin, end, prods2new_ents, 1);
        offset = end;
      }
    }
    old_lows2new_lows = old_ents2new_ents;
  }
  *mesh = new_mesh;
}

void amr_refine(Mesh* mesh, Bytes elems_are_marked) {
  OMEGA_H_CHECK(mesh->family() == OMEGA_H_HYPERCUBE);
  mark_amr(mesh, elems_are_marked);
  mesh->set_parting(OMEGA_H_GHOSTED);
  amr_refine_ghosted(mesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  amr_refine_elem_based(mesh);
}

}
