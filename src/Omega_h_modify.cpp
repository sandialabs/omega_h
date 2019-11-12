#include "Omega_h_modify.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_atomics.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_globals.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_owners.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_unmap_mesh.hpp"

namespace Omega_h {

static void modify_conn(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs prod_verts2verts, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents, LOs old_lows2new_lows) {
  OMEGA_H_TIME_FUNCTION;
  auto low_dim = ent_dim - 1;
  auto down_degree = element_degree(old_mesh->family(), ent_dim, low_dim);
  auto old_ents2old_lows = old_mesh->ask_down(ent_dim, low_dim);
  auto old_ent_lows2old_lows = old_ents2old_lows.ab2b;
  auto same_ent_lows2old_lows =
      unmap(same_ents2old_ents, old_ent_lows2old_lows, down_degree);
  auto same_ent_lows2new_lows =
      compound_maps(same_ent_lows2old_lows, old_lows2new_lows);
  auto prods2new_lows = Adj();
  if (low_dim > VERT) {
    auto new_low_verts2new_verts = new_mesh->ask_verts_of(low_dim);
    auto new_verts2new_lows = new_mesh->ask_up(VERT, low_dim);
    prods2new_lows = reflect_down(prod_verts2verts, new_low_verts2new_verts,
        new_verts2new_lows, old_mesh->family(), ent_dim, low_dim);
  } else {
    prods2new_lows = Adj(prod_verts2verts);
  }
  auto prod_lows2new_lows = prods2new_lows.ab2b;
  auto nsame_ents = same_ents2old_ents.size();
  auto nprods = prods2new_ents.size();
  auto nnew_ents = nsame_ents + nprods;
  auto ent_name = dimensional_singular_name(ent_dim);
  auto lows_name = dimensional_plural_name(low_dim);
  auto ent_lows2lows_name =
      std::string(ent_name) + " " + lows_name + " to " + lows_name;
  Write<LO> new_ent_lows2new_lows(nnew_ents * down_degree, ent_lows2lows_name);
  auto new_ent_low_codes = Write<I8>();
  map_into(
      prod_lows2new_lows, prods2new_ents, new_ent_lows2new_lows, down_degree);
  map_into(same_ent_lows2new_lows, same_ents2new_ents, new_ent_lows2new_lows,
      down_degree);
  if (low_dim > VERT) {
    auto codes_name = std::string(ent_name) + " " + lows_name + " codes";
    new_ent_low_codes = Write<I8>(nnew_ents * down_degree, codes_name);
    auto old_ent_low_codes = old_ents2old_lows.codes;
    OMEGA_H_CHECK(same_ents2old_ents.size() == same_ents2new_ents.size());
    auto same_ent_low_codes =
        read(unmap(same_ents2old_ents, old_ent_low_codes, down_degree));
    map_into(
        same_ent_low_codes, same_ents2new_ents, new_ent_low_codes, down_degree);
    auto prod_low_codes = prods2new_lows.codes;
    map_into(prod_low_codes, prods2new_ents, new_ent_low_codes, down_degree);
  }
  auto new_ents2new_lows =
      Adj(LOs(new_ent_lows2new_lows), Read<I8>(new_ent_low_codes));
  new_mesh->set_ents(ent_dim, new_ents2new_lows);
}

/* this is the case of AMR-style refinement where
   some entities being split (mods) are actually on
   the partition boundary.
   in this case they need to communicate on owners
   for the child entities.
   the rule is: a child copy is an owner if its parent
   copy is an owner */
static Remotes get_prod_owners_shared(Mesh* mesh, Few<LOs, 4> mods2mds,
    Few<LOs, 4> mods2prods, LOs prods2new_ents) {
  auto nprods = prods2new_ents.size();
  Write<I32> prod_own_ranks(nprods);
  Write<LO> prod_own_idxs(nprods);
  for (Int mod_dim = 0; mod_dim < 4; ++mod_dim) {
    if (!mods2prods[mod_dim].exists()) continue;
    auto prod_begin = mods2prods[mod_dim].first();
    auto prod_end = mods2prods[mod_dim].last();
    auto nmods = mods2prods[mod_dim].size() - 1;
    auto nmod_prods = prod_end - prod_begin;
    /* assuming the number of products per split entity is constant! */
    LO nprods_per_mod;
    if (nmods == 0) nprods_per_mod = 0;
    else nprods_per_mod = divide_no_remainder(nmod_prods, nmods);
    nprods_per_mod = mesh->comm()->allreduce(nprods_per_mod, OMEGA_H_MAX);
    auto md_ranks = mesh->ask_owners(mod_dim).ranks;
    auto mod_ranks = read(unmap(mods2mds[mod_dim], md_ranks, 1));
    auto mod_prod_idxs = unmap_range(prod_begin, prod_end, prods2new_ents, 1);
    mod_prod_idxs = mesh->sync_subset_array(
        mod_dim, mod_prod_idxs, mods2mds[mod_dim], -1, nprods_per_mod);
    map_into_range(mod_prod_idxs, prod_begin, prod_end, prod_own_idxs, 1);
    expand_into(mod_ranks, mods2prods[mod_dim], prod_own_ranks, 1);
  }
  return {prod_own_ranks, prod_own_idxs};
}

/* set the owners of the mesh after an adaptive rebuild pass.
   the entities that stay the same retain the same conceptual
   ownership, just updated by unmap_owners() to reflect new indices.
   all entities produced by this MPI rank exist only on this
   MPI rank, so they are their own owners */
static void modify_owners(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    Few<LOs, 4> mods2mds, Few<LOs, 4> mods2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents, LOs old_ents2new_ents,
    bool mods_can_be_shared) {
  OMEGA_H_TIME_FUNCTION;
  auto same_owners =
      unmap_owners(old_mesh, ent_dim, same_ents2old_ents, old_ents2new_ents);
  auto same_own_ranks = same_owners.ranks;
  auto same_own_idxs = same_owners.idxs;
  auto nprods = prods2new_ents.size();
  Remotes prod_owners;
  if (mods_can_be_shared) {
    prod_owners =
        get_prod_owners_shared(old_mesh, mods2mds, mods2prods, prods2new_ents);
  } else {
    prod_owners.ranks = Read<I32>(nprods, new_mesh->comm()->rank());
    prod_owners.idxs = prods2new_ents;
  }
  auto nsame_ents = same_ents2old_ents.size();
  auto nnew_ents = nsame_ents + nprods;
  Write<I32> new_own_ranks(nnew_ents);
  Write<LO> new_own_idxs(nnew_ents);
  map_into(same_own_ranks, same_ents2new_ents, new_own_ranks, 1);
  map_into(same_own_idxs, same_ents2new_ents, new_own_idxs, 1);
  map_into(prod_owners.ranks, prods2new_ents, new_own_ranks, 1);
  map_into(prod_owners.idxs, prods2new_ents, new_own_idxs, 1);
  auto new_owners = Remotes(Read<I32>(new_own_ranks), LOs(new_own_idxs));
  new_mesh->set_owners(ent_dim, new_owners);
}

static LOs collect_same(
    Mesh* mesh, Int ent_dim, Few<Bytes, 4> mds_are_mods, bool keep_mods) {
  if (keep_mods || (!mds_are_mods[ent_dim].exists())) {
    return LOs(mesh->nents(ent_dim), 0, 1);
  }
  auto ents_not_modified = invert_marks(mds_are_mods[ent_dim]);
  return collect_marked(ents_not_modified);
}

/* To determine the ordering (both global and local) of
   new entities of a certain dimension (ent_dim), each modified entity will
   be associated with a "representative entity" of dimension (ent_dim)
   in the old mesh.
   Conceptually, we will pretend that that representative entity
   "becomes" the entities "produced by" the modified entity.
   This function computes the mapping from modified entities to representatives.

   In some cases, two or more modified entities may have the same
   representative. This is the reason for the atomic_add in get_rep_counts(),
   the exch_reduce(SUM) later on, and the need to choose an ordering amongst
   modified entities that share a representative.
 */
LOs get_mods2reps(Mesh* mesh, Int ent_dim, Int mod_dim, LOs mods2mds) {
  auto nmods = mods2mds.size();
  LOs mods2reps;
  if (mod_dim == ent_dim) {
    mods2reps = mods2mds;
  } else if (ent_dim > mod_dim) {
    /* if the modified entity has lower dimension than the products,
       the first upward adjacent entity will represent
       this modified entity during the updating of global (and local) numbers.
       local ordering should be sorted by old globals
       already, so this is actually the adjacent entity
       with the lowest old (global and local) number */
    auto mds2ents = mesh->ask_up(mod_dim, ent_dim);
    Write<LO> mods2reps_w(nmods);
    auto setup_reps = OMEGA_H_LAMBDA(LO mod) {
      mods2reps_w[mod] = mds2ents.ab2b[mds2ents.a2ab[mods2mds[mod]]];
    };
    parallel_for(nmods, setup_reps, "get_mods2reps(up)");
    mods2reps = mods2reps_w;
  } else {
    /* if the modified entity has higher dimension than the products,
       the first downward adjacent entity will be the representative. */
    auto mds2ents = mesh->ask_down(mod_dim, ent_dim).ab2b;
    auto deg = element_degree(mesh->family(), mod_dim, ent_dim);
    Write<LO> mods2reps_w(nmods);
    auto setup_reps = OMEGA_H_LAMBDA(LO mod) {
      mods2reps_w[mod] = mds2ents[mods2mds[mod] * deg];
    };
    parallel_for(nmods, setup_reps, "get_mods2reps(down)");
    mods2reps = mods2reps_w;
  }
  return mods2reps;
}

/* return an array which maps an old entity (e)
   of dimension (ent_dim) to the number of new
   entities that (e) represents.
   This number is 1 for same entities, 0 for entities
   being removed, and
   N for representative of one or more modified entities, where N is
   the sum of the numbers of entities of dimension (ent_dim)
   which are produced by each modified entity, plus 1
   in the case where the representative is not going away.
   If (count_modified) is true, modified entities who represent themselves count
   themselves. If (count_non_owned) is false, non-owned entities do not count
   themselves.
 */
LOs get_rep_counts(Mesh* mesh, Int ent_dim, Few<LOs, 4> mods2mds,
    Few<LOs, 4> mods2reps, Few<LOs, 4> mods2nprods, LOs same_ents2ents,
    bool count_non_owned) {
  auto nents = mesh->nents(ent_dim);
  Write<LO> rep_counts(nents, 0);
  auto nsame_ents = same_ents2ents.size();
  auto ents_are_owned = mesh->owned(ent_dim);
  OMEGA_H_CHECK(ents_are_owned.size() == nents);
  auto mark_same = OMEGA_H_LAMBDA(LO same_ent) {
    auto ent = same_ents2ents[same_ent];
    if (count_non_owned || ents_are_owned[ent]) rep_counts[ent] = 1;
  };
  parallel_for(nsame_ents, mark_same, "get_rep_counts(same)");
  for (Int mod_dim = 0; mod_dim <= mesh->dim(); ++mod_dim) {
    if (!mods2reps[mod_dim].exists()) continue;
    auto mods2reps_dim = mods2reps[mod_dim];
    auto mods2nprods_dim = mods2nprods[mod_dim];
    auto mods2mds_dim = mods2mds[mod_dim];
    auto mds_are_owned = mesh->owned(mod_dim);
    auto nmods = mods2reps_dim.size();
    OMEGA_H_CHECK(nmods == mods2nprods_dim.size());
    auto mark_reps = OMEGA_H_LAMBDA(LO mod) {
      auto rep = mods2reps_dim[mod];
      auto nmod_prods = mods2nprods_dim[mod];
      /* TODO: the non-owned check can be skipped entirely for the AMR case */
      if (count_non_owned || mds_are_owned[mods2mds_dim[mod]]) {
        /* TODO: the atomic_add can be skipped for mod_dim <= ent_dim */
        atomic_add(&rep_counts[rep], nmod_prods);
      }
    };
    parallel_for(nmods, mark_reps, "get_rep_counts(modified)");
  }
  return rep_counts;
}

/* One more thing we need to for the
   case when multiple modified entities share a representative.
   That representative represents the sum of all adjacent modifications
   and has a single global number from which the globals
   of all represented product entities will be derived.
   We just need to determine an order in which to
   product globals from the representative.
   The most intuitive order is the upward adjacent ordering
   from representatives to products (guaranteed to be sorted by globals).
   The ordering is stored as a tag on edges.
   Each representative entity assigns a value to each
   modified entity that it represents.
   The global version of this only works in OMEGA_H_GHOSTED mode,
   so this function is called to save that ordering while in OMEGA_H_GHOSTED
   mode, and that global ordering is used later by assign_new_numbering(), which
   runs in OMEGA_H_ELEM_BASED mode.
   This function is also called to determine the local version
   of this ordering, i.e. only for modified entities on the same
   MPI rank in OMEGA_H_ELEM_BASED mode.
   The local version is used to deterime new *local* numbers
   for products, and is also the reason
   why this function doesn't check for ghosting.
 */

LOs get_rep2md_order_adapt(
    Mesh* mesh, Int key_dim, Int rep_dim, Bytes kds_are_keys) {
  Few<LOs, 4> mods2mds;
  Few<LOs, 4> mods2nprods;
  mods2mds[key_dim] = collect_marked(kds_are_keys);
  mods2nprods[key_dim] = LOs(mods2mds[key_dim].size(), 1);
  Few<bool, 4> mods_have_prods = {false, false, false, false};
  mods_have_prods[key_dim] = true;
  auto all_orders =
      get_rep2md_order(mesh, rep_dim, mods2mds, mods2nprods, mods_have_prods);
  OMEGA_H_CHECK(all_orders[key_dim].exists());
  return all_orders[key_dim];
}

Few<LOs, 4> get_rep2md_order(Mesh* mesh, Int rep_dim, Few<LOs, 4> mods2mds,
    Few<LOs, 4> mods2nprods, Few<bool, 4> mods_have_prods) {
  auto nreps = mesh->nents(rep_dim);
  auto elem_dim = mesh->dim();
  Few<Write<LO>, 4> rep2md_order_w;
  Few<Adj, 4> reps2mds;
  Few<LOs, 4> mds2mods;
  bool do_anything = false;
  for (Int mod_dim = rep_dim + 1; mod_dim <= elem_dim; ++mod_dim) {
    if (!mods_have_prods[mod_dim])
      continue;
    else
      do_anything = true;
    rep2md_order_w[mod_dim] = Write<LO>(mesh->nents(mod_dim), -1);
    reps2mds[mod_dim] = mesh->ask_up(rep_dim, mod_dim);
    mds2mods[mod_dim] =
        invert_injective_map(mods2mds[mod_dim], mesh->nents(mod_dim));
  }
  Few<LOs, 4> out;
  if (!do_anything) return out;
  if (mods_have_prods[rep_dim]) {
    mds2mods[rep_dim] =
        invert_injective_map(mods2mds[rep_dim], mesh->nents(rep_dim));
  }
  auto f = OMEGA_H_LAMBDA(LO rep) {
    LO offset = 0;
    if (mods_have_prods[rep_dim]) {
      auto mod = mds2mods[rep_dim][rep];
      if (mod >= 0) {
        offset += mods2nprods[rep_dim][mod];
      }
    }
    for (Int mod_dim = rep_dim + 1; mod_dim <= elem_dim; ++mod_dim) {
      if (!mods_have_prods[mod_dim]) continue;
      for (auto rep_md = reps2mds[mod_dim].a2ab[rep];
           rep_md < reps2mds[mod_dim].a2ab[rep + 1]; ++rep_md) {
        auto md = reps2mds[mod_dim].ab2b[rep_md];
        auto mod = mds2mods[mod_dim][md];
        if (mod >= 0) {
          auto code = reps2mds[mod_dim].codes[rep_md];
          auto dir = code_which_down(code);
          if (dir == 0) {
            rep2md_order_w[mod_dim][md] = offset;
            offset += mods2nprods[mod_dim][mod];
          }
        }
      }
    }
  };
  parallel_for(nreps, f, "get_rep2md_order");
  for (Int mod_dim = rep_dim + 1; mod_dim <= elem_dim; ++mod_dim) {
    if (!mods_have_prods[mod_dim]) continue;
    out[mod_dim] = LOs(rep2md_order_w[mod_dim]);
  }
  return out;
}

/* Assigns a new numbering to all entities in the new mesh,
   given the results of the scan of annotated old entities
   (old_ents2new_offsets). Because the construction of the new mesh depends on
   these numbers, we are still dealing with things separated into entities that
   stay the same, and newly produced entities. This function is mainly
   responsible for numbering newly produced entities based on the number that
   their representative entity got from the scan. */
template <typename T>
void assign_new_numbering(Read<T> old_ents2new_numbers,
    LOs same_ents2old_ents, Few<LOs, 4> mods2mds, Few<LOs, 4> mods2reps,
    Few<LOs, 4> mods2prods, Few<LOs, 4> rep2md_order,
    Read<T>* p_same_ents2new_numbers, Read<T>* p_prods2new_numbers,
    bool keep_mods) {
  *p_same_ents2new_numbers = unmap(same_ents2old_ents, old_ents2new_numbers, 1);
  LO nprods = 0;
  for (Int mod_dim = 0; mod_dim < 4; ++mod_dim) {
    if (mods2prods[mod_dim].exists()) nprods = mods2prods[mod_dim].last();
  }
  Write<T> prods2new_offsets_w(nprods);
  for (Int mod_dim = 0; mod_dim < 4; ++mod_dim) {
    if (!mods2prods[mod_dim].exists()) continue;
    auto mods2new_offsets = unmap(mods2reps[mod_dim], old_ents2new_numbers, 1);
    auto nmods = mods2reps[mod_dim].size();
    OMEGA_H_CHECK(nmods == mods2prods[mod_dim].size() - 1);
    Int rep_self_count;
    if (keep_mods) {
      /* currently we only keep modified entities in the case of AMR refinement
       */
      rep_self_count = 1;
    } else {
      /* if modified entities aren't being kept, I assume we're doing simplex
         adaptation, in which case representatives counting themselves should
         only happen if we're splitting simplex edges. This is a HACK */
      rep_self_count = (rep2md_order[mod_dim].exists() ? 1 : 0);
    }
    auto mods2prods_dim = mods2prods[mod_dim];
    auto mods2mds_dim = mods2mds[mod_dim];
    if (rep2md_order[mod_dim].exists()) {
      /* this is the upward adjacency case */
      OMEGA_H_CHECK(mods2mds[mod_dim].exists());
      auto rep2md_order_dim = rep2md_order[mod_dim];
      auto write_prod_offsets = OMEGA_H_LAMBDA(LO mod) {
        auto md = mods2mds_dim[mod];
        auto md_order = rep2md_order_dim[md];
        OMEGA_H_CHECK(md_order >= 0);
        auto offset = mods2new_offsets[mod] + md_order + rep_self_count;
        for (auto prod = mods2prods_dim[mod]; prod < mods2prods_dim[mod + 1];
             ++prod) {
          prods2new_offsets_w[prod] = offset++;
        }
      };
      parallel_for(nmods, std::move(write_prod_offsets));
    } else {
      auto write_prod_offsets = OMEGA_H_LAMBDA(LO mod) {
        auto offset = mods2new_offsets[mod] + rep_self_count;
        for (auto prod = mods2prods_dim[mod]; prod < mods2prods_dim[mod + 1];
             ++prod) {
          prods2new_offsets_w[prod] = offset++;
        }
      };
      parallel_for(nmods, std::move(write_prod_offsets));
    }
  }
  *p_prods2new_numbers = prods2new_offsets_w;
}

static void modify_globals(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    bool keep_mods, Few<LOs, 4> mods2mds, Few<LOs, 4> mods2prods,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents,
    Few<LOs, 4> mods2reps, LOs global_rep_counts) {
  OMEGA_H_TIME_FUNCTION;
  auto const nsame_ents = same_ents2old_ents.size();
  OMEGA_H_CHECK(nsame_ents == same_ents2new_ents.size());
  auto const nprods = prods2new_ents.size();
  auto const old_globals = old_mesh->globals(ent_dim);
  auto const comm = old_mesh->comm();
  auto const old_ents2lins = copies_to_linear_owners(comm, old_globals);
  auto const lins2old_ents = old_ents2lins.invert();
  auto const nlins = lins2old_ents.nroots();
  auto const lin_rep_counts =
      old_ents2lins.exch_reduce(global_rep_counts, 1, OMEGA_H_SUM);
  OMEGA_H_CHECK(lin_rep_counts.size() == nlins);
  auto const lin_globals = rescan_globals(old_mesh, lin_rep_counts);
  auto const old_ents2new_globals =
      lins2old_ents.exch(Read<GO>(lin_globals), 1);
  Few<LOs, 4> global_rep2md_order;
  for (Int mod_dim = 0; mod_dim <= old_mesh->dim(); ++mod_dim) {
    if ((mod_dim > ent_dim) && mods2prods[mod_dim].exists()) {
      auto const name =
          std::string("rep_") + hypercube_singular_name(ent_dim) + "2md_order";
      global_rep2md_order[mod_dim] = old_mesh->get_array<LO>(mod_dim, name);
    }
  }
  Read<GO> same_ents2new_globals;
  Read<GO> prods2new_globals;
  assign_new_numbering(old_ents2new_globals, same_ents2old_ents, mods2mds,
      mods2reps, mods2prods, global_rep2md_order, &same_ents2new_globals,
      &prods2new_globals, keep_mods);
  auto nnew_ents = new_mesh->nents(ent_dim);
  OMEGA_H_CHECK(nnew_ents == nsame_ents + nprods);
  Write<GO> new_globals(nnew_ents);
  map_into(same_ents2new_globals, same_ents2new_ents, new_globals, 1);
  map_into(prods2new_globals, prods2new_ents, new_globals, 1);
  new_mesh->add_tag(ent_dim, "global", 1, Read<GO>(new_globals));
}

void modify_ents_adapt(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prod_verts2verts, LOs old_lows2new_lows,
    LOs* p_prods2new_ents, LOs* p_same_ents2old_ents, LOs* p_same_ents2new_ents,
    LOs* p_old_ents2new_ents) {
  OMEGA_H_TIME_FUNCTION;
  Few<LOs, 4> mods2mds;
  Few<Bytes, 4> mds_are_mods;
  Few<LOs, 4> mods2prods;
  mods2mds[key_dim] = keys2kds;
  mds_are_mods[key_dim] = mark_image(keys2kds, old_mesh->nents(key_dim));
  for (Int mod_dim = key_dim + 1; mod_dim <= old_mesh->dim(); ++mod_dim) {
    mds_are_mods[mod_dim] =
        mark_up(old_mesh, key_dim, mod_dim, mds_are_mods[key_dim]);
    mods2mds[mod_dim] = collect_marked(mds_are_mods[mod_dim]);
  }
  mods2prods[key_dim] = keys2prods;
  modify_ents(old_mesh, new_mesh, ent_dim, mods2mds, mds_are_mods, mods2prods,
      prod_verts2verts, old_lows2new_lows,
      /*keep_mods*/ false, /*mods_can_be_shared*/ false, p_prods2new_ents,
      p_same_ents2old_ents, p_same_ents2new_ents, p_old_ents2new_ents);
}

void modify_ents(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    Few<LOs, 4> mods2mds, Few<Bytes, 4> mds_are_mods, Few<LOs, 4> mods2prods,
    LOs prod_verts2verts, LOs old_lows2new_lows, bool keep_mods,
    bool mods_can_be_shared, LOs* p_prods2new_ents, LOs* p_same_ents2old_ents,
    LOs* p_same_ents2new_ents, LOs* p_old_ents2new_ents) {
  OMEGA_H_TIME_FUNCTION;
  *p_same_ents2old_ents =
      collect_same(old_mesh, ent_dim, mds_are_mods, keep_mods);
  Few<LOs, 4> mods2nprods;
  Few<LOs, 4> mods2reps;
  for (Int mod_dim = 0; mod_dim <= old_mesh->dim(); ++mod_dim) {
    if (!mods2prods[mod_dim].exists()) continue;
    OMEGA_H_CHECK(mods2prods[mod_dim].size() == mods2mds[mod_dim].size() + 1);
    mods2reps[mod_dim] =
        get_mods2reps(old_mesh, ent_dim, mod_dim, mods2mds[mod_dim]);
    mods2nprods[mod_dim] = get_degrees(mods2prods[mod_dim]);
  }
  auto local_rep_counts = get_rep_counts(old_mesh, ent_dim, mods2mds, mods2reps,
      mods2nprods, *p_same_ents2old_ents, /*count_non_owned*/ true);
  auto local_offsets = offset_scan(local_rep_counts);
  auto nnew_ents = local_offsets.last();
  Few<bool, 4> mods_have_prods;
  for (Int mod_dim = 0; mod_dim < 4; ++mod_dim)
    mods_have_prods[mod_dim] = mods2prods[mod_dim].exists();
  auto local_rep2md_order = get_rep2md_order(
      old_mesh, ent_dim, mods2mds, mods2nprods, mods_have_prods);
  assign_new_numbering(local_offsets, *p_same_ents2old_ents, mods2mds,
      mods2reps, mods2prods, local_rep2md_order, p_same_ents2new_ents,
      p_prods2new_ents, keep_mods);
  auto nold_ents = old_mesh->nents(ent_dim);
  *p_old_ents2new_ents =
      map_onto(*p_same_ents2new_ents, *p_same_ents2old_ents, nold_ents, -1, 1);
  if (ent_dim == VERT) {
    new_mesh->set_verts(nnew_ents);
  } else {
    modify_conn(old_mesh, new_mesh, ent_dim, prod_verts2verts,
        *p_prods2new_ents, *p_same_ents2old_ents, *p_same_ents2new_ents,
        old_lows2new_lows);
  }
  if (old_mesh->comm()->size() > 1) {
    modify_owners(old_mesh, new_mesh, ent_dim, mods2mds, mods2prods,
        *p_prods2new_ents, *p_same_ents2old_ents, *p_same_ents2new_ents,
        *p_old_ents2new_ents, mods_can_be_shared);
  }
  auto global_rep_counts = get_rep_counts(old_mesh, ent_dim, mods2mds,
      mods2reps, mods2nprods, *p_same_ents2old_ents, /*count_non_owned*/ false);
  modify_globals(old_mesh, new_mesh, ent_dim, keep_mods, mods2mds, mods2prods,
      *p_prods2new_ents, *p_same_ents2old_ents, *p_same_ents2new_ents,
      mods2reps, global_rep_counts);
}

void set_owners_by_indset(
    Mesh* mesh, Int key_dim, LOs keys2kds, Graph kds2elems) {
  if (mesh->comm()->size() == 1) return;
  auto kd_owners = mesh->ask_owners(key_dim);
  auto nkeys = keys2kds.size();
  auto elem_dim = mesh->dim();
  auto kds2kd_elems = kds2elems.a2ab;
  auto kd_elems2elems = kds2elems.ab2b;
  auto elem_owners = mesh->ask_owners(elem_dim);
  auto elems2owners = mesh->ask_dist(elem_dim);
  auto new_elem_ranks = deep_copy(elem_owners.ranks);
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto kd = keys2kds[key];
    auto kd_rank = kd_owners.ranks[kd];
    for (auto kd_elem = kds2kd_elems[kd]; kd_elem < kds2kd_elems[kd + 1];
         ++kd_elem) {
      auto elem = kd_elems2elems[kd_elem];
      new_elem_ranks[elem] = kd_rank;
    }
  };
  parallel_for(nkeys, f, "set_owners_by_indset");
  elem_owners = update_ownership(elems2owners, new_elem_ranks);
  mesh->set_owners(elem_dim, elem_owners);
}

}  // end namespace Omega_h
