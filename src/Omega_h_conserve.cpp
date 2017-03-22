#include "transfer_conserve.hpp"

namespace Omega_h {

void transfer_conserve_refine(Mesh* old_mesh, XferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_conserve(old_mesh, opts, dim, tagbase)) {
      transfer_inherit_refine(old_mesh, new_mesh, keys2edges, old_mesh->dim(),
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase);
    }
  }
}

template <Int dim>
static void transfer_conserve_dim(Mesh* old_mesh, Mesh* new_mesh,
    TagBase const* tagbase, Graph keys2old_mat_elems, Graph keys2new_mat_elems,
    Write<Real> new_data_w) {
  auto ncomps = tagbase->ncomps();
  auto old_tag = to<Real>(tagbase);
  auto old_data = old_tag->array();
  auto new_data = new_mesh->get_array<Real>(dim, tagbase->name());
  auto old_sizes = old_mesh->ask_sizes();
  auto new_sizes = new_mesh->ask_sizes();
}

static void transfer_conserve_tag(Mesh* old_mesh, Mesh* new_mesh,
    std::vector<std::pair<Graph, Graph>> mats_keys2elems,
    LOs same_ents2old_ents, LOs same_ents2new_ents, TagBase const* tagbase) {
  auto dim = old_mesh->dim();
  auto nnew_elems = new_mesh->nelems();
  auto ncomps = tagbase->ncomps();
  auto new_data_w = Write<Real>(nnew_elems * ncomps);
  for (auto pair : mats_keys2elems) {
    auto keys2old_mat_elems = pair.first;
    auto keys2new_mat_elems = pair.second;
    if (dim == 3)
      transfer_conserve_dim<3>(old_mesh, new_mesh, tagbase, keys2old_mat_elems,
          keys2new_mat_elems, new_data_w);
    if (dim == 2)
      transfer_conserve_dim<2>(old_mesh, new_mesh, tagbase, keys2old_mat_elems,
          keys2new_mat_elems, new_data_w);
  }
  transfer_common2(old_mesh, new_mesh, dim, same_ents2old_ents,
      same_ents2new_ents, tagbase, new_data_w);
}

void transfer_conserve(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
    Int key_dim, LOs keys2kds, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  auto dim = new_mesh->dim();
  if (!should_conserve_any(old_mesh, opts)) return;
  auto kds2old_elems = old_mesh->ask_up(key_dim, dim);
  auto keys2old_elems = unmap_graph(keys2kds, kds2old_elems);
  auto keys2new_elems = Graph(keys2prods, prods2new_ents);
  std::vector<std::pair<Graph, Graph>> mats_keys2elems;
  if (old_mesh->has_tag(dim, "class_id")) {
    auto old_class_ids = old_mesh->get_array<I32>(dim, "class_id");
    auto new_class_ids = new_mesh->get_array<I32>(dim, "class_id");
    mats_keys2elems = separate_cavities(
        keys2old_elems, old_class_ids, keys2new_elems, new_class_ids);
  } else {
    mats_keys2elems.push_back({keys2old_elems, keys2new_elems});
  }
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_conserve(old_mesh, opts, dim, tagbase)) {
      transfer_conserve_tag(old_mesh, new_mesh, mats_keys2elems,
          same_ents2old_ents, same_ents2new_ents, tagbase);
    }
  }
}

void fix_momentum_velocity_verts(
    Mesh* mesh, Int class_dim, I32 class_id, Int comp) {
  for (Int ent_dim = VERT; ent_dim <= class_dim; ++ent_dim) {
    auto ent_marks = mark_class_closure(mesh, ent_dim, class_dim, class_id);
    auto comp_marks = multiply_each_by(I8(1 << comp), ent_marks);
    if (mesh->has_tag(ent_dim, "momentum_velocity_fixed")) {
      auto old_marks = mesh->get_array<I8>(ent_dim, "momentum_velocity_fixed");
      auto new_marks = bit_or_each(old_marks, comp_marks);
      mesh->set_tag(ent_dim, "momentum_velocity_fixed", new_marks);
    } else {
      mesh->add_tag(ent_dim, "momentum_velocity_fixed", 1, comp_marks);
    }
  }
  for (Int ent_dim = class_dim + 1; ent_dim <= mesh->dim(); ++ent_dim) {
    if (!mesh->has_tag(ent_dim, "momentum_velocity_fixed")) {
      mesh->add_tag(ent_dim, "momentum_velocity_fixed", 1,
          Read<I8>(mesh->nents(ent_dim), I8(0)));
    }
  }
}

bool has_fixed_momentum_velocity(Mesh* mesh, XferOpts const& opts) {
  return has_momentum_velocity(mesh, opts) &&
         mesh->has_tag(VERT, "momentum_velocity_fixed");
}

Read<I8> filter_coarsen_momentum_velocity(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  auto comps_are_fixed = mesh->get_array<I8>(VERT, "momentum_velocity_fixed");
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto ncands = cands2edges.size();
  auto out = Write<I8>(ncands);
  auto f = LAMBDA(LO cand) {
    auto e = cands2edges[cand];
    auto code = cand_codes[cand];
    for (Int eev_col = 0; eev_col < 2; ++eev_col) {
      if (!collapses(code, eev_col)) continue;
      auto v_col = ev2v[e * 2 + eev_col];
      I8 fixed_comps = 7;
      for (auto ve = v2e.a2ab[v_col]; ve < v2e.a2ab[v_col + 1]; ++ve) {
        auto e2 = v2e.ab2b[ve];
        auto e2_code = v2e.codes[ve];
        auto eev_in = code_which_down(e2_code);
        auto eev_out = 1 - eev_in;
        auto ov = ev2v[e2 * 2 + eev_out];
        fixed_comps &= comps_are_fixed[ov];
      }
      if (fixed_comps) code = dont_collapse(code, eev_col);
    }
    out[cand] = code;
  };
  parallel_for(ncands, f);
  return mesh->sync_subset_array(
      EDGE, Read<I8>(out), cands2edges, I8(DONT_COLLAPSE), 1);
}

Read<I8> filter_swap_momentum_velocity(Mesh* mesh, LOs cands2edges) {
  auto vert_comps_are_fixed =
      mesh->get_array<I8>(VERT, "momentum_velocity_fixed");
  auto edges2elems = mesh->ask_up(EDGE, mesh->dim());
  auto cands2elems = unmap_graph(cands2edges, edges2elems);
  auto elem_verts2vert = mesh->ask_elem_verts();
  auto nverts_per_elem = simplex_degrees[mesh->dim()][VERT];
  constexpr Int max_verts_per_cavity = (AvgDegree<3, 0, 3>::value + 1) * 2;
  auto ncands = cands2elems.nnodes();
  auto keep_w = Write<I8>(ncands);
  auto f = LAMBDA(LO cand) {
    Few<LO, max_verts_per_cavity> cavity_verts;
    Int ncavity_verts = 0;
    for (auto ce = cands2elems.a2ab[cand]; ce < cands2elems.a2ab[cand + 1];
         ++ce) {
      auto elem = cands2elems.ab2b[ce];
      for (auto ev = elem * nverts_per_elem; ev < (elem + 1) * nverts_per_elem;
           ++ev) {
        auto vert = elem_verts2vert[ev];
        add_unique(cavity_verts, ncavity_verts, vert);
      }
    }
    Int cavity_bits = (1 << 3) - 1;
    for (Int i = 0; i < ncavity_verts; ++i) {
      auto vert = cavity_verts[i];
      auto vert_bits = Int(vert_comps_are_fixed[vert]);
      cavity_bits &= vert_bits;
    }
    keep_w[cand] = (cavity_bits == 0);
  };
  parallel_for(ncands, f);
  auto keep = Read<I8>(keep_w);
  return mesh->sync_subset_array(EDGE, keep, cands2edges, I8(0), 1);
}

static Reals get_cavity_momenta(
    Mesh* mesh, Graph keys2elems, Reals vert_velocities, Reals elem_masses) {
  auto dim = mesh->dim();
  auto cavity_elem_masses = unmap(keys2elems.ab2b, elem_masses, 1);
  auto cavity_elem_velocities =
      average_field(mesh, dim, keys2elems.ab2b, dim, vert_velocities);
  auto cavity_elem_momenta =
      multiply_each(cavity_elem_velocities, cavity_elem_masses);
  return fan_reduce(keys2elems.a2ab, cavity_elem_momenta, dim, OMEGA_H_SUM);
}

static Read<I8> get_comps_are_fixed(Mesh* mesh) {
  if (mesh->has_tag(VERT, "momentum_velocity_fixed")) {
    return mesh->get_array<I8>(VERT, "momentum_velocity_fixed");
  } else {
    return Read<I8>(mesh->nverts(), I8(0));
  }
}

template <Int dim>
void momentum_velocity_part1_dim(Mesh* donor_mesh, XferOpts const& opts,
    Mesh* target_mesh, Int key_dim, LOs keys2kds, LOs keys2prods,
    LOs prods2new_elems) {
  if (!has_momentum_velocity(donor_mesh, opts)) return;
  auto nkeys = keys2kds.size();
  auto keys2target_elems = Graph(keys2prods, prods2new_elems);
  auto kds2donor_elems = donor_mesh->ask_up(key_dim, dim);
  auto keys2donor_elems = unmap_graph(keys2kds, kds2donor_elems);
  auto target_elems2verts = target_mesh->ask_elem_verts();
  auto comps_are_fixed = get_comps_are_fixed(target_mesh);
  for (Int tag_i = 0; tag_i < donor_mesh->ntags(VERT); ++tag_i) {
    auto tagbase = donor_mesh->get_tag(VERT, tag_i);
    if (!is_momentum_velocity(donor_mesh, opts, VERT, tagbase)) continue;
    CHECK(tagbase->ncomps() == dim);
    auto tag = to<Real>(tagbase);
    auto donor_vert_velocities = tag->array();
    auto mass_name = opts.momentum_map.find(tag->name())->second;
    auto donor_elem_masses = donor_mesh->get_array<Real>(dim, mass_name);
    auto target_elem_masses = target_mesh->get_array<Real>(dim, mass_name);
    auto donor_cavity_momenta = get_cavity_momenta(
        donor_mesh, keys2donor_elems, donor_vert_velocities, donor_elem_masses);
    auto target_vert_velocities =
        target_mesh->get_array<Real>(VERT, tag->name());
    auto target_cavity_momenta = get_cavity_momenta(target_mesh,
        keys2target_elems, target_vert_velocities, target_elem_masses);
    auto cavity_momentum_losses =
        subtract_each(donor_cavity_momenta, target_cavity_momenta);
    auto corrections_w =
        Write<Real>(target_mesh->nelems() * (dim + 1) * dim, 0.0);
    auto f = LAMBDA(LO key) {
      auto lost_momentum = get_vector<dim>(cavity_momentum_losses, key);
      Few<Real, dim> free_momentum;
      Few<Int, dim> nfree_uses;
      for (Int comp = 0; comp < dim; ++comp) {
        free_momentum[comp] = 0;
        nfree_uses[comp] = 0;
      }
      for (auto ke = keys2target_elems.a2ab[key];
           ke < keys2target_elems.a2ab[key + 1]; ++ke) {
        auto e = keys2target_elems.ab2b[ke];
        auto e_mass = target_elem_masses[e];
        for (auto eev = 0; eev < dim + 1; ++eev) {
          auto v = target_elems2verts[e * (dim + 1) + eev];
          auto v_vel = get_vector<dim>(target_vert_velocities, v);
          auto eev_momentum = v_vel * (e_mass / (dim + 1));
          auto code = Int(comps_are_fixed[v]);
          for (Int comp = 0; comp < dim; ++comp) {
            if (!((1 << comp) & code)) {
              ++nfree_uses[comp];
              free_momentum[comp] += fabs(eev_momentum[comp]);
            }
          }
        }
      }
      for (auto ke = keys2target_elems.a2ab[key];
           ke < keys2target_elems.a2ab[key + 1]; ++ke) {
        auto e = keys2target_elems.ab2b[ke];
        auto e_mass = target_elem_masses[e];
        for (auto eev = 0; eev < dim + 1; ++eev) {
          auto v = target_elems2verts[e * (dim + 1) + eev];
          auto v_vel = get_vector<dim>(target_vert_velocities, v);
          auto eev_momentum = v_vel * (e_mass / (dim + 1));
          auto code = Int(comps_are_fixed[v]);
          for (Int comp = 0; comp < dim; ++comp) {
            if (!((1 << comp) & code)) {
              Real comp_corr;
              if (free_momentum[comp] > EPSILON) {
                comp_corr = lost_momentum[comp] *
                            (fabs(eev_momentum[comp]) / free_momentum[comp]);
              } else {
                comp_corr =
                    lost_momentum[comp] * (1.0 / Real(nfree_uses[comp]));
              }
              corrections_w[(e * (dim + 1) + eev) * dim + comp] = comp_corr;
            }
          }
        }
      }
    };
    parallel_for(nkeys, f);
    auto name = tag->name() + "_correction";
    target_mesh->add_tag(dim, name, (dim + 1) * dim, Reals(corrections_w));
  }
}

void do_momentum_velocity_part1(Mesh* donor_mesh, XferOpts const& opts,
    Mesh* target_mesh, Int key_dim, LOs keys2kds, LOs keys2prods,
    LOs prods2new_elems) {
  if (donor_mesh->dim() == 3) {
    momentum_velocity_part1_dim<3>(donor_mesh, opts, target_mesh, key_dim,
        keys2kds, keys2prods, prods2new_elems);
  } else if (donor_mesh->dim() == 2) {
    momentum_velocity_part1_dim<2>(donor_mesh, opts, target_mesh, key_dim,
        keys2kds, keys2prods, prods2new_elems);
  } else {
    NORETURN();
  }
}

static Reals get_vertex_masses(Mesh* mesh, Reals elems2mass) {
  CHECK(mesh->owners_have_all_upward(VERT));
  auto dim = mesh->dim();
  auto verts2elems = mesh->ask_up(VERT, dim);
  auto verts2mass = graph_reduce(verts2elems, elems2mass, 1, OMEGA_H_SUM);
  verts2mass = multiply_each_by(1.0 / Real(dim + 1), verts2mass);
  return mesh->sync_array(VERT, verts2mass, 1);
}

template <Int dim>
void momentum_velocity_part2_dim(Mesh* mesh, XferOpts const& opts) {
  if (!has_momentum_velocity(mesh, opts)) return;
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto v2e = mesh->ask_up(VERT, dim);
  auto comps_are_fixed = get_comps_are_fixed(mesh);
  for (Int tag_i = 0; tag_i < mesh->ntags(VERT); ++tag_i) {
    auto tagbase = mesh->get_tag(VERT, tag_i);
    if (!is_momentum_velocity(mesh, opts, VERT, tagbase)) continue;
    CHECK(tagbase->ncomps() == dim);
    auto tag = to<Real>(tagbase);
    auto corr_name = tag->name() + "_correction";
    auto mass_name = opts.momentum_map.find(tag->name())->second;
    auto elem_masses = mesh->get_array<Real>(dim, mass_name);
    auto vert_masses = get_vertex_masses(mesh, elem_masses);
    auto elem_corrections = mesh->get_array<Real>(dim, corr_name);
    auto vert_corrections_w = Write<Real>(mesh->nverts() * dim, 0.0);
    auto f = LAMBDA(LO v) {
      Vector<dim> correction = zero_vector<dim>();
      for (auto ve = v2e.a2ab[v]; ve < v2e.a2ab[v + 1]; ++ve) {
        auto e = v2e.ab2b[ve];
        auto code = v2e.codes[ve];
        auto eev = code_which_down(code);
        for (Int comp = 0; comp < dim; ++comp) {
          correction[comp] +=
              elem_corrections[(e * (dim + 1) + eev) * dim + comp];
        }
      }
      set_vector(vert_corrections_w, v, correction);
    };
    parallel_for(mesh->nverts(), f);
    auto vert_corrections = Reals(vert_corrections_w);
    vert_corrections = mesh->sync_array(VERT, vert_corrections, dim);
    auto old_velocities = tag->array();
    auto velocity_corrections = divide_each(vert_corrections, vert_masses);
    auto new_velocities = add_each(old_velocities, velocity_corrections);
    mesh->set_tag(VERT, tag->name(), new_velocities);
    mesh->remove_tag(dim, corr_name);
  }
}

void do_momentum_velocity_part2(Mesh* mesh, XferOpts const& opts) {
  if (mesh->dim() == 3)
    momentum_velocity_part2_dim<3>(mesh, opts);
  else if (mesh->dim() == 2)
    momentum_velocity_part2_dim<2>(mesh, opts);
  else
    NORETURN();
}

}  // end namespace Omega_h
