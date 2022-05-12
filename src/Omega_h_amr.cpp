#include <Omega_h_amr.hpp>
#include <Omega_h_amr_topology.hpp>
#include <Omega_h_amr_transfer.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_globals.hpp>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_int_scan.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_modify.hpp>
#include <Omega_h_unmap_mesh.hpp>

namespace Omega_h {

namespace amr {

void remove_non_leaf_uses(Mesh* mesh) {
  Bytes ent_persists[4];
  auto elem_dim = mesh->dim();
  ent_persists[elem_dim] = mesh->ask_leaves(elem_dim);
  for (Int ent_dim = 0; ent_dim < elem_dim; ++ent_dim) {
    ent_persists[ent_dim] = mark_down(
        mesh, elem_dim, ent_dim, ent_persists[elem_dim]);
  }
  LOs new_ents2old_ents[4];
  for (Int ent_dim = 0; ent_dim <= elem_dim; ++ent_dim) {
    new_ents2old_ents[ent_dim] = collect_marked(ent_persists[ent_dim]);
  }
  unmap_mesh(mesh, new_ents2old_ents);
}

static OMEGA_H_DEVICE Byte should_elem_be_refined(LO elem, Adj elems2bridges,
    Adj bridges2elems, Bytes is_interior, Bytes is_bridge_leaf,
    Int nbridges_per_elem, Children children, Bytes elems_are_marked,
    Write<Byte> one_level_mark) {
  Byte mark = 0;
  for (Int b = 0; b < nbridges_per_elem; ++b) {
    auto bridge = elems2bridges.ab2b[elem * nbridges_per_elem + b];
    if (!is_interior[bridge]) continue;
    if (is_bridge_leaf[bridge]) continue;
    auto bridge_child_begin = children.a2ab[bridge];
    auto bridge_child_end = children.a2ab[bridge + 1];
    for (auto c = bridge_child_begin; c < bridge_child_end; ++c) {
      auto child = children.ab2b[c];
      auto child_adj_elem_begin = bridges2elems.a2ab[child];
      auto child_adj_elem_end = bridges2elems.a2ab[child + 1];
      OMEGA_H_CHECK((child_adj_elem_end - child_adj_elem_begin) == 1);
      auto child_adj_elem = bridges2elems.ab2b[child_adj_elem_begin];
      if (elems_are_marked[child_adj_elem]) mark = 1;
      if (one_level_mark[child_adj_elem]) mark = 1;
    }
  }
  return mark;
}

Bytes enforce_2to1_refine(Mesh* mesh, Int bridge_dim, Bytes elems_are_marked) {
  auto elem_dim = mesh->dim();
  OMEGA_H_CHECK(bridge_dim > 0);
  OMEGA_H_CHECK(bridge_dim < elem_dim);
  auto is_elem_leaf = mesh->ask_leaves(elem_dim);
  auto is_bridge_leaf = mesh->ask_leaves(bridge_dim);
  auto elems2bridges = mesh->ask_down(elem_dim, bridge_dim);
  auto bridges2elems = mesh->ask_up(bridge_dim, elem_dim);
  auto nbridges_per_elem = Omega_h::hypercube_degree(elem_dim, bridge_dim);
  auto is_interior = Omega_h::mark_by_class_dim(mesh, bridge_dim, elem_dim);
  auto children = mesh->ask_children(bridge_dim, bridge_dim);
  Write<Byte> one_level_mark(mesh->nelems());
  auto f = OMEGA_H_LAMBDA(LO elem) {
    if (!is_elem_leaf[elem]) {
      one_level_mark[elem] = 0;
    } else if (elems_are_marked[elem]) {
      one_level_mark[elem] = 1;
    } else {
      one_level_mark[elem] = should_elem_be_refined(elem, elems2bridges,
          bridges2elems, is_interior, is_bridge_leaf, nbridges_per_elem,
          children, elems_are_marked, one_level_mark);
    }
  };
  Omega_h::parallel_for(mesh->nelems(), f, "enforce_one_level");
  return one_level_mark;
}

static void refine_ghosted(Mesh* mesh) {
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
    auto rep2md_orders = get_rep2md_order(
        mesh, prod_dim, mods2mds, mods2nprods, mods_have_prods);
    auto name =
        std::string("rep_") + hypercube_singular_name(prod_dim) + "2md_order";
    for (Int mod_dim = prod_dim + 1; mod_dim <= mesh->dim(); ++mod_dim) {
      mesh->add_tag(mod_dim, name, 1, rep2md_orders[mod_dim]);
    }
  }
}

static void refine_elem_based(Mesh* mesh, TransferOpts xfer_opts) {
  auto prod_counts = amr::count_refined(mesh);
  Few<Bytes, 4> mds_are_mods;
  Few<LOs, 4> mods2mds;
  Few<LOs, 4> mds2mods;
  for (Int mod_dim = 0; mod_dim <= mesh->dim(); ++mod_dim) {
    mds_are_mods[mod_dim] = mesh->get_array<Byte>(mod_dim, "refine");
    mods2mds[mod_dim] = collect_marked(mds_are_mods[mod_dim]);
    mds2mods[mod_dim] =
        invert_injective_map(mods2mds[mod_dim], mesh->nents(mod_dim));
  }
  auto new_mesh = mesh->copy_meta();
  Few<LOs, 4> mods2midverts;
  Few<LOs, 4> prods2new_ents;
  Few<LOs, 4> same_ents2old_ents;
  Few<LOs, 4> same_ents2new_ents;
  Few<LOs, 4> old_ents2new_ents;
  LOs old_lows2new_lows;
  for (Int prod_dim = 0; prod_dim <= mesh->dim(); ++prod_dim) {
    LOs prods2verts;
    if (prod_dim != VERT) {
      prods2verts =
          amr::get_refined_topology(mesh, prod_dim, prod_counts[prod_dim],
              mods2mds, mds2mods, mods2midverts, old_ents2new_ents[0]);
    }
    Few<LOs, 4> mods2prods;
    {
      LO offset = 0;
      for (Int mod_dim = max2(Int(EDGE), prod_dim); mod_dim <= mesh->dim();
           ++mod_dim) {
        auto nprods_per_mod = hypercube_split_degree(mod_dim, prod_dim);
        auto nmods_of_dim = mods2mds[mod_dim].size();
        mods2prods[mod_dim] =
            LOs(mods2mds[mod_dim].size() + 1, offset, nprods_per_mod);
        offset += nprods_per_mod * nmods_of_dim;
      }
    }
    modify_ents(mesh, &new_mesh, prod_dim, mods2mds, mds_are_mods, mods2prods,
        prods2verts, old_lows2new_lows, /*keep_mods*/ true,
        /*mods_can_be_shared*/ true, &(prods2new_ents[prod_dim]),
        &(same_ents2old_ents[prod_dim]), &(same_ents2new_ents[prod_dim]),
        &(old_ents2new_ents[prod_dim]));
    if (prod_dim == VERT) {
      mods2midverts[VERT] =
          unmap(mods2mds[VERT], old_ents2new_ents[prod_dim], 1);
      LO offset = 0;
      for (Int mod_dim = EDGE; mod_dim <= mesh->dim(); ++mod_dim) {
        OMEGA_H_CHECK(hypercube_split_degree(mod_dim, prod_dim) == 1);
        auto nmods_of_dim = mods2mds[mod_dim].size();
        auto begin = offset;
        auto end = begin + nmods_of_dim;
        mods2midverts[mod_dim] =
            unmap_range(begin, end, prods2new_ents[prod_dim], 1);
        offset = end;
      }
      amr::transfer_linear_interp(mesh, &new_mesh, mods2mds, mods2midverts,
          same_ents2old_ents[prod_dim], same_ents2new_ents[prod_dim],
          xfer_opts);
    }
    amr::transfer_levels(mesh, &new_mesh, prod_dim, mods2mds,
        prods2new_ents[prod_dim], same_ents2old_ents[prod_dim],
        same_ents2new_ents[prod_dim]);
    amr::transfer_leaves(mesh, &new_mesh, prod_dim, mods2mds,
        prods2new_ents[prod_dim], same_ents2old_ents[prod_dim],
        same_ents2new_ents[prod_dim], old_ents2new_ents[prod_dim]);
    old_lows2new_lows = old_ents2new_ents[prod_dim];
  }
  amr::transfer_parents(mesh, &new_mesh, mods2mds, prods2new_ents,
      same_ents2old_ents, same_ents2new_ents, old_ents2new_ents);
  amr::transfer_inherit(mesh, &new_mesh, prods2new_ents, same_ents2old_ents,
      same_ents2new_ents, xfer_opts);
  *mesh = new_mesh;
}

void refine(Mesh* mesh, Bytes elems_are_marked, TransferOpts xfer_opts) {
  OMEGA_H_CHECK(mesh->family() == OMEGA_H_HYPERCUBE);
  amr::mark_refined(mesh, elems_are_marked);
  mesh->set_parting(OMEGA_H_GHOSTED);
  amr::refine_ghosted(mesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  amr::refine_elem_based(mesh, xfer_opts);
}

void derefine(Mesh* mesh, Bytes elems_are_marked, TransferOpts xfer_opts) {
  OMEGA_H_CHECK(mesh->family() == OMEGA_H_HYPERCUBE);
  amr::tag_derefined(mesh, elems_are_marked);
  LOs new_ents2old_ents[4];
  GOs new_globals[4];
  for (int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    auto does_ent_persist = mesh->get_array<Byte>(ent_dim, "persists");
    new_ents2old_ents[ent_dim] = collect_marked(does_ent_persist);
    mesh->remove_tag(ent_dim, "persists");
    auto old_ents2new_globals = rescan_globals(mesh, does_ent_persist);
    new_globals[ent_dim] =
        unmap(new_ents2old_ents[ent_dim], old_ents2new_globals, 1);
  }
  unmap_mesh(mesh, new_ents2old_ents);
  for (int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    mesh->set_tag(ent_dim, "global", new_globals[ent_dim]);
  }
  (void)xfer_opts;
}

}  // namespace amr

}  // namespace Omega_h
