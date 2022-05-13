#include "Omega_h_migrate.hpp"

#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_owners.hpp"

namespace Omega_h {

Remotes form_down_use_owners(Mesh* mesh, Int high_dim, Int low_dim) {
  auto uses2lows = mesh->ask_down(high_dim, low_dim).ab2b;
  auto lows2owners = mesh->ask_owners(low_dim);
  return unmap(uses2lows, lows2owners);
}

Dist get_old_owners2uniq_uses(Dist uses2old_owners) {
  auto old_owners2uses = uses2old_owners.invert();
  auto nold_owners = old_owners2uses.nroots();
  auto nserv_uses = old_owners2uses.nitems();
  auto serv_uses2ranks = old_owners2uses.items2ranks();
  auto old_owners2serv_uses = old_owners2uses.roots2items();
  Write<I8> keep_w(nserv_uses, 1);
  Write<LO> degrees_w(nold_owners);
  auto f = OMEGA_H_LAMBDA(LO old_owner) {
    LO degree = 0;
    auto begin = old_owners2serv_uses[old_owner];
    auto end = old_owners2serv_uses[old_owner + 1];
    for (LO serv_use = begin; serv_use < end; ++serv_use) {
      if (!keep_w[serv_use]) continue;
      ++degree;
      auto rank = serv_uses2ranks[serv_use];
      for (LO serv_use2 = serv_use + 1; serv_use2 < end; ++serv_use2) {
        if (serv_uses2ranks[serv_use2] == rank) {
          keep_w[serv_use2] = 0;
        }
      }
    }
    degrees_w[old_owner] = degree;
  };
  parallel_for(nold_owners, f, "get_old_owners2uniq_uses");
  auto degrees = LOs(degrees_w);
  auto keep = Read<I8>(keep_w);
  auto uniq_serv_uses2serv_uses = collect_marked(keep);
  auto uniq_serv_uses2ranks =
      unmap(uniq_serv_uses2serv_uses, serv_uses2ranks, 1);
  auto old_owners2uniq_serv_uses = offset_scan(degrees);
  Dist old_owners2uniq_uses;
  old_owners2uniq_uses.set_parent_comm(uses2old_owners.parent_comm());
  old_owners2uniq_uses.set_dest_ranks(uniq_serv_uses2ranks);
  old_owners2uniq_uses.set_roots2items(old_owners2uniq_serv_uses);
  return old_owners2uniq_uses;
}

Dist get_new_copies2old_owners(Dist uses2old_owners, GOs old_owner_globals) {
  auto old_owners2uniq_uses = get_old_owners2uniq_uses(uses2old_owners);
  auto old_owners2items = old_owners2uniq_uses.roots2items();
  old_owners2uniq_uses.set_dest_globals(old_owner_globals);
  auto uniq_uses2old_owners = old_owners2uniq_uses.invert();
  return uniq_uses2old_owners;
}

LOs form_new_conn(Dist new_ents2old_owners, Dist old_owners2new_uses) {
  auto nnew_ents = new_ents2old_owners.nitems();
  auto serv_ents2new_idxs = new_ents2old_owners.exch(LOs(nnew_ents, 0, 1), 1);
  auto old_owners2new_ents = new_ents2old_owners.invert();
  auto serv_ents2ranks = old_owners2new_ents.items2ranks();
  auto serv_uses2ranks = old_owners2new_uses.items2ranks();
  auto nserv_uses = old_owners2new_uses.nitems();
  auto old_owners2serv_uses = old_owners2new_uses.roots2items();
  auto old_owners2serv_ents = old_owners2new_ents.roots2items();
  auto nold_owners = old_owners2new_ents.nroots();
  Write<LO> serv_uses2new_idxs(nserv_uses);
  auto f = OMEGA_H_LAMBDA(LO old_owner) {
    auto ebegin = old_owners2serv_ents[old_owner];
    auto eend = old_owners2serv_ents[old_owner + 1];
    auto ubegin = old_owners2serv_uses[old_owner];
    auto uend = old_owners2serv_uses[old_owner + 1];
    for (auto serv_use = ubegin; serv_use < uend; ++serv_use) {
      auto rank = serv_uses2ranks[serv_use];
      LO idx = -1;
      for (auto e = ebegin; e < eend; ++e) {
        if (serv_ents2ranks[e] == rank) {
          idx = serv_ents2new_idxs[e];
          break;
        }
      }
      serv_uses2new_idxs[serv_use] = idx;
    }
  };
  parallel_for(nold_owners, f, "form_new_conn");
  auto serv_uses2new_uses = old_owners2new_uses;
  serv_uses2new_uses.set_roots2items(LOs());
  return serv_uses2new_uses.exch(LOs(serv_uses2new_idxs), 1);
}

void push_down(Mesh* old_mesh, Int ent_dim, Int low_dim,
    Dist old_owners2new_ents, Adj& new_ents2new_lows,
    Dist& old_low_owners2new_lows) {
  OMEGA_H_TIME_FUNCTION;
  auto nlows_per_high = element_degree(old_mesh->family(), ent_dim, low_dim);
  auto old_use_owners = form_down_use_owners(old_mesh, ent_dim, low_dim);
  auto new_use_owners =
      old_owners2new_ents.exch(old_use_owners, nlows_per_high);
  Dist low_uses2old_owners(
      old_mesh->comm(), new_use_owners, old_mesh->nents(low_dim));
  auto old_owner_globals = old_mesh->globals(low_dim);
  auto new_lows2old_owners =
      get_new_copies2old_owners(low_uses2old_owners, old_owner_globals);
  old_low_owners2new_lows = new_lows2old_owners.invert();
  auto old_low_owners2new_uses = low_uses2old_owners.invert();
  auto new_conn = form_new_conn(new_lows2old_owners, old_low_owners2new_uses);
  new_ents2new_lows.ab2b = new_conn;
  auto old_codes = old_mesh->ask_down(ent_dim, low_dim).codes;
  if (old_codes.exists()) {
    auto new_codes = old_owners2new_ents.exch(old_codes, nlows_per_high);
    new_ents2new_lows.codes = new_codes;
  }
}

void push_tags(Mesh *old_mesh, Mesh* new_mesh, Int ent_dim,
    Dist old_owners2new_ents) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(old_owners2new_ents.nroots() == old_mesh->nents(ent_dim));
  ScopedChangeRCFieldsToMesh rc_to_mesh(*old_mesh);
  for (Int i = 0; i < old_mesh->ntags(ent_dim); ++i) {
    auto tag = old_mesh->get_tag(ent_dim, i);
    auto const& name = tag->name();
    const auto ncomps = tag->ncomps();
    const auto class_ids = tag->class_ids();
    detail::apply_to_omega_h_types(tag->type(), [&](auto t){
      using T = decltype(t);
      auto array = old_mesh->get_array<T>(ent_dim,name);
      array = old_owners2new_ents.exch(array, ncomps);

      if(is_rc_tag(name) && rc_to_mesh.did_conversion() ) {
        new_mesh->set_rc_from_mesh_array(ent_dim,ncomps,class_ids,name,array);
      }
      else {
      // FIXME missing class_ids for rc tag
        new_mesh->add_tag<T>(ent_dim, name, ncomps, array, true);
      }

    });
  }
}

void push_ents(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    Dist new_ents2old_owners, Dist old_owners2new_ents, Omega_h_Parting mode) {
  push_tags(old_mesh, new_mesh, ent_dim, old_owners2new_ents);
  Read<I32> own_ranks;
  /* if we are ghosting, each entity should remain owned by the
   * same rank that owned it before ghosting, as this is the only
   * mechanism we have to identify the ghost layers.
   * if we are doing a vertex-based partitioning, at least the
   * vertices ought to retain their original owners, for similar
   * reasons.
   */
  if ((mode == OMEGA_H_GHOSTED) ||
      ((mode == OMEGA_H_VERT_BASED) && (ent_dim == VERT))) {
    auto old_own_ranks = old_mesh->ask_owners(ent_dim).ranks;
    own_ranks = old_owners2new_ents.exch(old_own_ranks, 1);
  }
  auto owners = update_ownership(new_ents2old_owners, own_ranks);
  new_mesh->set_owners(ent_dim, owners);
}

static void print_migrate_stats(CommPtr comm, Dist new_elems2old_owners) {
  auto msgs2ranks = new_elems2old_owners.msgs2ranks();
  auto msgs2content = new_elems2old_owners.msgs2content();
  auto msg_content_size = get_degrees(msgs2content);
  auto msg_content_size_h = HostRead<LO>(msg_content_size);
  auto msgs2ranks_h = HostRead<I32>(msgs2ranks);
  LO nintern = 0;
  for (LO msg = 0; msg < msgs2ranks_h.size(); ++msg) {
    if (msgs2ranks_h[msg] == comm->rank()) {
      nintern = msg_content_size_h[msg];
    }
  }
  auto npulled = msgs2content.last();
  auto nextern = npulled - nintern;
  auto total_pulled = comm->allreduce(GO(npulled), OMEGA_H_SUM);
  auto total_extern = comm->allreduce(GO(nextern), OMEGA_H_SUM);
  if (comm->rank() == 0) {
    std::cout << "migration pulling (" << total_extern << " remote) / ("
              << total_pulled << " total) elements\n";
  }
}

void migrate_matches(Mesh* mesh, Mesh* new_mesh, Int const d,
    Dist const* old_owners2new_ents) {
  auto size = mesh->comm()->size();
  auto rank = mesh->comm()->rank();
  auto r2i = old_owners2new_ents->roots2items();//roots2items
  auto i2dR = old_owners2new_ents->items2ranks();//items2destRanks
  auto i2dI = old_owners2new_ents->items2dest_idxs();//items2destIds
  auto nents = mesh->nents(d);
  std::vector<int> dest_r;//ranks
  std::vector<int> dest_i;//idxs
  HostWrite<LO> new_r2i(nents+1, 0, 0, "froots2fitems");
  HostWrite<LO> old_isMatch(nents, -1, 0, "old_isMatch");
  //roots2items stores offsets; +1 for offsets. TODO:it is possible to
  //make this proportional to O(max [froot IDs of all(i.e.
  //missing+present) rroots])
  if (nents) {//if mesh exists on a process
    auto matches = mesh->get_matches(d);
    HostRead<LO> leaf_idxs(matches.leaf_idxs);
    HostRead<LO> root_idxs(matches.root_idxs);
    HostRead<LO> root_ranks(matches.root_ranks);
    HostRead<LO> owners_idxs((mesh->ask_owners(d)).idxs);
    HostRead<LO> h_r2i(r2i);
    HostRead<LO> h_i2dR(i2dR);
    HostRead<LO> h_i2dI(i2dI);
    auto max_leaf = leaf_idxs.last();
    auto max_leaf_rootID = root_idxs.last();
    auto max_leaf_rootRank = root_ranks.last();
    if (max_leaf_rootRank == rank) {
      OMEGA_H_CHECK(max_leaf_rootID != max_leaf);
    //check leaf and root are disinct
    }
    
    LO ent = 0;
    for (LO i = 0; i < matches.leaf_idxs.size(); ++i) {
      auto leaf = leaf_idxs[i];
      auto root = root_idxs[i];
      auto root_rank = root_ranks[i];

      /*for missing old ents to avoid broken graphs and hanging nodes*/
      auto leaf_next = leaf;
      while (ent < leaf_next) {
        leaf = ent;
        LO root_owner = 0;
        LO n_items = 0;
        root_owner = owners_idxs[leaf];//no matching hence root=leaf=ent
        auto L_R_item_begin = h_r2i[root_owner];
        auto L_R_item_end = h_r2i[root_owner+1];
        auto L_R = L_R_item_end - L_R_item_begin;
        for (int item = L_R_item_begin; item < L_R_item_end; ++item) {
          auto L_rank = h_i2dR[item];
          auto L_idx = h_i2dI[item];
          dest_r.push_back(L_rank);
          dest_i.push_back(L_idx);
        }
        n_items += L_R;
        new_r2i[root_owner+1] = n_items;
        ++ent;
      }
      ent = leaf_next+1;
      leaf = leaf_next;
      /**/

      if (((rank == root_rank) && (leaf < root)) || (rank != root_rank)) {
      //distributed i/p: if ((leaf < root)||(rank < root_rank))
        LO root_owner = 0;
        //for distributed mesh input as the roots will
        //be present on different parts hence we cannot get the info about the
        //root's new copies. TODO: create dist from
        //oldLeafOwners-2-oldRootOwners to exch info about their new copies
        root_owner = owners_idxs[root];
        auto L_R_item_begin = h_r2i[root_owner];
        auto L_R_item_end = h_r2i[root_owner+1];
        auto L_R = L_R_item_end - L_R_item_begin;
        auto leaf_owner = owners_idxs[leaf];
        auto L_L_item_begin = h_r2i[leaf_owner];
        auto L_L_item_end = h_r2i[leaf_owner+1];
        auto L_L = L_L_item_end - L_L_item_begin;
        auto L = L_R + L_L;
        const LO n_items = L;
        new_r2i[leaf_owner+1] = n_items;//leaf_owner cause for picking new
        //owners the anchor should be one having lower old GID
        old_isMatch[leaf_owner] = 1;
        for (int item = L_R_item_begin; item < L_R_item_end; ++item) {
          auto L_rank = h_i2dR[item];
          auto L_idx = h_i2dI[item];
          dest_r.push_back(L_rank);
          dest_i.push_back(L_idx);
        }
        for (int item = L_L_item_begin; item < L_L_item_end; ++item) {
          auto L_rank = h_i2dR[item];
          auto L_idx = h_i2dI[item];
          dest_r.push_back(L_rank);
          dest_i.push_back(L_idx);
        }
      }
    }
    //create offsets from nitems
    LO n_items = 0;
    for (LO j = 1; j < new_r2i.size(); ++j) {
      n_items += new_r2i[j];
      new_r2i[j] = n_items;
    }
  } // end if mesh has ents
  HostWrite<I32> host_dest_r(dest_r.size());
  HostWrite<LO> host_dest_i(dest_r.size());
  for (unsigned int i = 0; i < dest_r.size(); ++i) {
    host_dest_r[i] = dest_r[static_cast<std::size_t>(i)];
    host_dest_i[i] = dest_i[static_cast<std::size_t>(i)];
  }

  HostWrite<LO> max_dest_id(size, -1, 0, "max_dest_id");
  if (nents) {
    for (LO leaf = 0; leaf < host_dest_i.size(); ++leaf) {
      I32 destR = host_dest_r[leaf];
      LO destId = host_dest_i[leaf];
      if (destId > max_dest_id[destR]) max_dest_id[destR] = destId;
    }
  }
  std::vector<int> mesh_dest_r;//destination ranks for mesh (all ranks)
  std::vector<int> max_dest_i;//destination ranks for mesh (all ranks)
  if (nents) {
    for (LO i = 0; i<size; ++i) {
      mesh_dest_r.push_back(i);
      max_dest_i.push_back(max_dest_id[i]);
    }
  }
  HostWrite<LO> mesh_dest_r_h(mesh_dest_r.size());
  HostWrite<LO> max_dest_i_h(max_dest_i.size());
  for (unsigned int i = 0; i < mesh_dest_r.size(); ++i) {
    mesh_dest_r_h[i] = mesh_dest_r[static_cast<std::size_t>(i)];
    max_dest_i_h[i] = max_dest_i[static_cast<std::size_t>(i)];
  }
  Dist old_rank2new_ranks;
  old_rank2new_ranks.set_parent_comm(mesh->comm());
  old_rank2new_ranks.set_dest_ranks(mesh_dest_r_h.write());

  auto rev_max_dest_ids =
    old_rank2new_ranks.exch(Read<LO>(max_dest_i_h.write()), 1);
  //old_rank2new_ranks.exch_reduce(Read<LO>(max_dest_i_h.write()), 1,
  //OMEGA_H_MAX);//TODO: for parallel input we want to get max
  //of dest ids recved from all old ranks
  LO rev_nroots;
  rev_nroots = rev_max_dest_ids.last()+1;

  Dist owners2new_leaves;
  owners2new_leaves.set_parent_comm(mesh->comm());
  owners2new_leaves.set_dest_ranks(host_dest_r.write());
  int wait=0; while(wait);
  owners2new_leaves.set_dest_idxs(host_dest_i.write(), rev_nroots);
  owners2new_leaves.set_roots2items(new_r2i.write());
  Read<I32> own_ranks;
  auto new_matchOwners = update_ownership(owners2new_leaves.invert(), own_ranks);
  auto new_isMatch = owners2new_leaves.exch(Read<LO>(old_isMatch.write()), 1);
  HostRead<LO> new_isMatch_h(new_isMatch);

  /*create c_r in new mesh*/
  std::vector<int> leaf_ids;
  std::vector<int> root_ids;
  std::vector<int> root_rks;
  HostRead<LO> new_matchOwner_idxs(new_matchOwners.idxs);
  HostRead<I32> new_matchOwner_ranks(new_matchOwners.ranks);

  for (LO i = 0; i < new_matchOwner_idxs.size(); ++i) {
    if (new_isMatch_h[i] == 1) {
      leaf_ids.push_back(i);
      root_ids.push_back(new_matchOwner_idxs[i]);
      root_rks.push_back(new_matchOwner_ranks[i]);
    }
  }
  HostWrite<LO> host_leaf_ids(leaf_ids.size());
  HostWrite<LO> host_root_ids(leaf_ids.size());
  HostWrite<LO> host_root_rks(leaf_ids.size());
  for (unsigned int i = 0; i < leaf_ids.size(); ++i) {
    host_leaf_ids[i] = leaf_ids[static_cast<std::size_t>(i)]; 
    host_root_ids[i] = root_ids[static_cast<std::size_t>(i)]; 
    host_root_rks[i] = root_rks[static_cast<std::size_t>(i)]; 
  }
  auto cr = c_Remotes(LOs(host_leaf_ids.write()),
                      LOs(host_root_ids.write()),
                      LOs(host_root_rks.write()));
  new_mesh->set_matches(d, cr);
}

void migrate_mesh(
    Mesh* mesh, Dist new_elems2old_owners, Omega_h_Parting mode, bool verbose) {
  OMEGA_H_TIME_FUNCTION;
  for (Int d = 0; d <= mesh->dim(); ++d) {
    OMEGA_H_CHECK(mesh->has_tag(d, "global"));
  }
  auto new_mesh = mesh->copy_meta();
  auto comm = mesh->comm();
  auto dim = mesh->dim();
  if (verbose) print_migrate_stats(comm, new_elems2old_owners);
  Dist new_ents2old_owners = new_elems2old_owners;
  auto old_owners2new_ents = new_ents2old_owners.invert();
  for (Int d = dim; d > VERT; --d) {
    Adj high2low;
    Dist old_low_owners2new_lows;
    push_down(
        mesh, d, d - 1, old_owners2new_ents, high2low, old_low_owners2new_lows);
    new_mesh.set_ents(d, high2low);
    new_ents2old_owners = old_owners2new_ents.invert();
    push_ents(
        mesh, &new_mesh, d, new_ents2old_owners, old_owners2new_ents, mode);

    if ((mesh->is_matched() > 0) && (d < dim)) {
      migrate_matches(mesh, &new_mesh, d, &old_owners2new_ents);
    }

    old_owners2new_ents = old_low_owners2new_lows;
  }

  auto new_verts2old_owners = old_owners2new_ents.invert();
  auto nnew_verts = new_verts2old_owners.nitems();
  new_mesh.set_verts(nnew_verts);
  push_ents(
      mesh, &new_mesh, VERT, new_verts2old_owners, old_owners2new_ents, mode);

  if (mesh->is_matched() > 0) {
    migrate_matches(mesh, &new_mesh, VERT, &old_owners2new_ents);
  }

  *mesh = new_mesh;
  for (Int d = 0; d <= mesh->dim(); ++d) {
    OMEGA_H_CHECK(mesh->has_tag(d, "global"));
  }

}

}  // end namespace Omega_h
