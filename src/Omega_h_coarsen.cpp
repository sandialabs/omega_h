#include "Omega_h_coarsen.hpp"

#include <iostream>

#include <cstdlib>
#include "Omega_h_array_ops.hpp"
#include "Omega_h_collapse.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_indset.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_modify.hpp"
#include "Omega_h_transfer.hpp"

namespace Omega_h {

static Read<I8> get_edge_codes(Mesh* mesh) {
  auto edge_cand_codes = mesh->get_array<I8>(EDGE, "collapse_code");
  mesh->remove_tag(EDGE, "collapse_code");
  return edge_cand_codes;
}

static void put_edge_codes(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  auto edge_codes =
      map_onto(cand_codes, cands2edges, mesh->nedges(), I8(DONT_COLLAPSE), 1);
  mesh->add_tag(EDGE, "collapse_code", 1, edge_codes);
}

static bool coarsen_element_based1(Mesh* mesh) {
  auto comm = mesh->comm();
  auto edge_cand_codes = get_edge_codes(mesh);
  auto edges_are_cands = each_neq_to(edge_cand_codes, I8(DONT_COLLAPSE));
  auto cands2edges = collect_marked(edges_are_cands);
  auto cand_codes = read(unmap(cands2edges, edge_cand_codes, 1));
  cand_codes = check_collapse_class(mesh, cands2edges, cand_codes);
  /* edge and endpoints classification check */
  if (get_max(comm, cand_codes) <= DONT_COLLAPSE) return false;
  put_edge_codes(mesh, cands2edges, cand_codes);
  return true;
}

static void filter_coarsen_candidates(
    LOs* cands2edges, Read<I8>* cand_codes, Reals* cand_quals = nullptr) {
  auto keep = each_neq_to(*cand_codes, I8(DONT_COLLAPSE));
  auto new2old = collect_marked(keep);
  *cands2edges = unmap(new2old, *cands2edges, 1);
  *cand_codes = unmap(new2old, *cand_codes, 1);
  if (cand_quals) *cand_quals = unmap(new2old, *cand_quals, 2);
}

enum OvershootLimit { DESIRED, ALLOWED };

enum Improve { DONT_IMPROVE, IMPROVE_LOCALLY };

static bool coarsen_ghosted(Mesh* mesh, AdaptOpts const& opts,
    OvershootLimit overshoot, Improve improve) {
  auto comm = mesh->comm();
  auto edge_cand_codes = get_edge_codes(mesh);
  auto edges_are_cands = each_neq_to(edge_cand_codes, I8(DONT_COLLAPSE));
  auto cands2edges = collect_marked(edges_are_cands);
  auto cand_edge_codes = read(unmap(cands2edges, edge_cand_codes, 1));
  /* surface exposure (classification) checks */
  cand_edge_codes = check_collapse_exposure(mesh, cands2edges, cand_edge_codes);
  filter_coarsen_candidates(&cands2edges, &cand_edge_codes);
  /* edge length overshoot check */
  auto max_length = (overshoot == DESIRED) ? opts.max_length_desired
                                           : opts.max_length_allowed;
  cand_edge_codes =
      prevent_coarsen_overshoot(mesh, max_length, cands2edges, cand_edge_codes);
  filter_coarsen_candidates(&cands2edges, &cand_edge_codes);
  if (comm->reduce_and(cands2edges.size() == 0)) return false;
  #ifndef _MSC_VER
  if (opts.should_prevent_coarsen_flip) {
    cand_edge_codes = prevent_coarsen_flip(mesh, cands2edges, cand_edge_codes);
    filter_coarsen_candidates(&cands2edges, &cand_edge_codes);
    if (comm->reduce_and(cands2edges.size() == 0)) return false;
  }
  #endif
  /* cavity quality checks */
  auto cand_edge_quals = coarsen_qualities(mesh, cands2edges, cand_edge_codes);
  cand_edge_codes = filter_coarsen_min_qual(
      cand_edge_codes, cand_edge_quals, opts.min_quality_allowed);
  if (improve == IMPROVE_LOCALLY) {
    cand_edge_codes = filter_coarsen_improve(
        mesh, cands2edges, cand_edge_codes, cand_edge_quals);
  }
  filter_coarsen_candidates(&cands2edges, &cand_edge_codes, &cand_edge_quals);
  /* finished cavity quality checks */
  if (comm->reduce_and(cands2edges.size() == 0)) return false;
  auto verts_are_cands = Read<I8>();
  auto vert_quals = Reals();
  auto vert_rails = Read<GO>();
  choose_rails(mesh, cands2edges, cand_edge_codes, cand_edge_quals,
      &verts_are_cands, &vert_quals, &vert_rails);
  auto verts_are_keys = find_indset(mesh, VERT, vert_quals, verts_are_cands);
  Graph verts2cav_elems;
  verts2cav_elems = mesh->ask_up(VERT, mesh->dim());
  mesh->add_tag(VERT, "key", 1, verts_are_keys);
  mesh->add_tag(VERT, "collapse_quality", 1, vert_quals);
  mesh->add_tag(VERT, "collapse_rail", 1, vert_rails);
  auto keys2verts = collect_marked(verts_are_keys);
  set_owners_by_indset(mesh, VERT, keys2verts, verts2cav_elems);
  return true;
}

static void coarsen_element_based2(Mesh* mesh, AdaptOpts const& opts) {
  auto comm = mesh->comm();
  auto verts_are_keys = mesh->get_array<I8>(VERT, "key");
  auto vert_quals = mesh->get_array<Real>(VERT, "collapse_quality");
  auto vert_rails = mesh->get_array<GO>(VERT, "collapse_rail");
  mesh->remove_tag(VERT, "collapse_rail");
  auto keys2verts = collect_marked(verts_are_keys);
  auto nkeys = keys2verts.size();
  if (opts.verbosity >= EACH_REBUILD) {
    auto ntotal_keys = comm->allreduce(GO(nkeys), OMEGA_H_SUM);
    if (comm->rank() == 0) {
      std::cout << "coarsening " << ntotal_keys << " vertices\n";
    }
  }
  auto rails2edges = LOs();
  auto rail_col_dirs = Read<I8>();
  find_rails(mesh, keys2verts, vert_rails, &rails2edges, &rail_col_dirs);
  auto dead_ents = mark_dead_ents(mesh, rails2edges, rail_col_dirs);
  auto keys2verts_onto = get_verts_onto(mesh, rails2edges, rail_col_dirs);
  auto new_mesh = mesh->copy_meta();
  auto old_verts2new_verts = LOs();
  auto old_lows2new_lows = LOs();
  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    auto keys2prods = LOs();
    auto prod_verts2verts = LOs();
    auto keys2doms = Adj();
    if (ent_dim == VERT) {
      keys2prods = LOs(nkeys + 1, 0);
    } else {
      keys2doms =
          find_coarsen_domains(mesh, keys2verts, ent_dim, dead_ents[ent_dim]);
      keys2prods = keys2doms.a2ab;
      prod_verts2verts = coarsen_topology(
          mesh, keys2verts_onto, ent_dim, keys2doms, old_verts2new_verts);
    }
    auto prods2new_ents = LOs();
    auto same_ents2old_ents = LOs();
    auto same_ents2new_ents = LOs();
    auto old_ents2new_ents = LOs();
    modify_ents_adapt(mesh, &new_mesh, ent_dim, VERT, keys2verts, keys2prods,
        prod_verts2verts, old_lows2new_lows, &prods2new_ents,
        &same_ents2old_ents, &same_ents2new_ents, &old_ents2new_ents);
    if (ent_dim == VERT) {
      old_verts2new_verts = old_ents2new_ents;
    }
    transfer_coarsen(mesh, opts.xfer_opts, &new_mesh, keys2verts, keys2doms,
        ent_dim, prods2new_ents, same_ents2old_ents, same_ents2new_ents);
    old_lows2new_lows = old_ents2new_ents;
  }
  *mesh = new_mesh;
}

static bool coarsen(Mesh* mesh, AdaptOpts const& opts, OvershootLimit overshoot,
    Improve improve) {
  begin_code("coarsen");
  auto ret = coarsen_element_based1(mesh);
  if (ret) {
    mesh->set_parting(OMEGA_H_GHOSTED);
    ret = coarsen_ghosted(mesh, opts, overshoot, improve);
  }
  if (ret) {
    mesh->set_parting(OMEGA_H_ELEM_BASED, false);
    coarsen_element_based2(mesh, opts);
  }
  end_code();
  return ret;
}

bool coarsen_verts(Mesh* mesh, AdaptOpts const& opts,
    Read<I8> vert_marks, OvershootLimit overshoot, Improve improve) {
  auto ev2v = mesh->ask_verts_of(EDGE);
  Write<I8> edge_codes_w(mesh->nedges(), DONT_COLLAPSE);
  auto f = OMEGA_H_LAMBDA(LO e) {
    I8 code = DONT_COLLAPSE;
    for (Int eev = 0; eev < 2; ++eev) {
      if (vert_marks[ev2v[e * 2 + eev]]) {
        code = do_collapse(code, eev);
      }
    }
    edge_codes_w[e] = code;
  };
  parallel_for(mesh->nedges(), f, "coarsen_verts(edge_codes)");
  mesh->add_tag(EDGE, "collapse_code", 1, Read<I8>(edge_codes_w));
  return coarsen(mesh, opts, overshoot, improve);
}

static bool coarsen_ents(Mesh* mesh, AdaptOpts const& opts, Int ent_dim,
    Read<I8> marks, OvershootLimit overshoot, Improve improve) {
  auto vert_marks = mark_down(mesh, ent_dim, VERT, marks);
  return coarsen_verts(mesh, opts, vert_marks, overshoot, improve);
}

bool coarsen_by_size(Mesh* mesh, AdaptOpts const& opts) {
  OMEGA_H_TIME_FUNCTION;
  auto comm = mesh->comm();
  auto lengths = mesh->ask_lengths();
  auto edge_is_cand = each_lt(lengths, opts.min_length_desired);
  auto ret = (get_max(comm, edge_is_cand) == 1);
  if (ret) {
    ret = coarsen_ents(mesh, opts, EDGE, edge_is_cand, DESIRED, DONT_IMPROVE);
  }
  return ret;
}

bool coarsen_slivers(Mesh* mesh, AdaptOpts const& opts) {
  OMEGA_H_TIME_FUNCTION;
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto comm = mesh->comm();
  auto elems_are_cands =
      mark_sliver_layers(mesh, opts.min_quality_desired, opts.nsliver_layers);
  OMEGA_H_CHECK(get_max(comm, elems_are_cands) == 1);
  auto ret = coarsen_ents(
      mesh, opts, mesh->dim(), elems_are_cands, ALLOWED, IMPROVE_LOCALLY);
  return ret;
}

}  // end namespace Omega_h
