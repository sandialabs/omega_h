#include "transfer_conserve.hpp"

#include "Omega_h_math.hpp"
#include "Omega_h_r3d.hpp"
#include "array.hpp"
#include "collapse.hpp"
#include "graph.hpp"
#include "loop.hpp"
#include "map.hpp"
#include "simplices.hpp"
#include "size.hpp"
#include "tag.hpp"
#include "transfer.hpp"

namespace Omega_h {

static void transfer_conserve_refine(Mesh* old_mesh, Mesh* new_mesh,
    LOs keys2edges, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents, std::string const& name) {
  auto prod_dim = old_mesh->dim();
  auto old_tag = old_mesh->get_tag<Real>(prod_dim, name);
  auto ncomps = old_tag->ncomps();
  auto nprods = keys2prods.last();
  auto prod_data = Write<Real>(nprods * ncomps);
  auto nkeys = keys2edges.size();
  /* transfer pairs */
  auto dom_dim = prod_dim;
  auto dom_data = old_mesh->get_array<Real>(dom_dim, name);
  auto edges2doms = old_mesh->ask_graph(EDGE, dom_dim);
  auto edges2edge_doms = edges2doms.a2ab;
  auto edge_doms2doms = edges2doms.ab2b;
  auto f = LAMBDA(LO key) {
    auto edge = keys2edges[key];
    auto prod = keys2prods[key];
    for (auto edge_dom = edges2edge_doms[edge];
         edge_dom < edges2edge_doms[edge + 1]; ++edge_dom) {
      auto dom = edge_doms2doms[edge_dom];
      for (Int pair = 0; pair < 2; ++pair) {
        for (Int comp = 0; comp < ncomps; ++comp) {
          prod_data[prod * ncomps + comp] = dom_data[dom * ncomps + comp] / 2.0;
        }
        ++prod;
      }
    }
  };
  parallel_for(nkeys, f);
  transfer_common(old_mesh, new_mesh, prod_dim, same_ents2old_ents,
      same_ents2new_ents, prods2new_ents, old_tag, Read<Real>(prod_data));
}

void transfer_conserve_refine(Mesh* old_mesh, Mesh* new_mesh, LOs keys2edges,
    LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (tagbase->xfer() == OMEGA_H_CONSERVE) {
      transfer_conserve_refine(old_mesh, new_mesh, keys2edges, keys2prods,
          prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase->name());
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
  auto old_ev2v = old_mesh->ask_elem_verts();
  auto old_coords = old_mesh->coords();
  auto new_ev2v = new_mesh->ask_elem_verts();
  auto new_coords = new_mesh->coords();
  auto measure = RealElementSizes(new_mesh);
  auto nkeys = keys2old_mat_elems.nnodes();
  enum { max_targets = AvgDegree<dim, 0, dim>::value * 2 };
  auto f = LAMBDA(LO key) {
    auto kte_begin = keys2new_mat_elems.a2ab[key];
    auto kte_end = keys2new_mat_elems.a2ab[key + 1];
    auto ntargets = kte_end - kte_begin;
    for (Int i = 0; i < ntargets; ++i) {
      auto target_elem = keys2new_mat_elems.ab2b[kte_begin + i];
      for (Int comp = 0; comp < ncomps; ++comp) {
        new_data_w[target_elem * ncomps + comp] = 0;
      }
    }
    for (auto kde = keys2old_mat_elems.a2ab[key];
         kde < keys2old_mat_elems.a2ab[key + 1]; ++kde) {
      auto donor_elem = keys2old_mat_elems.ab2b[kde];
      auto donor_verts = gather_verts<dim + 1>(old_ev2v, donor_elem);
      auto donor_points = gather_vectors<dim + 1, dim>(old_coords, donor_verts);
      Vector<max_targets> coeffs;
      Real total_size = 0.0;
      for (Int i = 0; i < ntargets; ++i) {
        auto target_elem = keys2new_mat_elems.ab2b[kte_begin + i];
        auto target_verts = gather_verts<dim + 1>(new_ev2v, target_elem);
        auto target_points =
            gather_vectors<dim + 1, dim>(new_coords, target_verts);
        auto intersection =
            r3d::intersect_simplices(target_points, donor_points);
        auto intersection_size = r3d::measure(intersection);
        coeffs[i] = intersection_size;
        total_size += intersection_size;
      }
      for (Int i = 0; i < ntargets; ++i) coeffs[i] /= total_size;
      for (Int i = 0; i < ntargets; ++i) {
        auto target_elem = keys2new_mat_elems.ab2b[kte_begin + i];
        for (Int comp = 0; comp < ncomps; ++comp) {
          new_data_w[target_elem * ncomps + comp] +=
              coeffs[i] * old_data[donor_elem * ncomps + comp];
        }
      }
    }
  };
  parallel_for(nkeys, f);
}

static void transfer_conserve_tag(Mesh* old_mesh, Mesh* new_mesh,
    std::map<Int, Graph> keys2old_elems_cat,
    std::map<Int, Graph> keys2new_elems_cat, LOs same_ents2old_ents,
    LOs same_ents2new_ents, TagBase const* tagbase) {
  auto dim = old_mesh->dim();
  auto nnew_elems = new_mesh->nelems();
  auto ncomps = tagbase->ncomps();
  auto new_data_w = Write<Real>(nnew_elems * ncomps);
  for (auto pair : keys2old_elems_cat) {
    auto mat = pair.first;
    auto keys2old_mat_elems = pair.second;
    CHECK(keys2new_elems_cat.count(mat));
    auto keys2new_mat_elems = keys2new_elems_cat[mat];
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

void transfer_conserve(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  auto dim = new_mesh->dim();
  if (!has_xfer(old_mesh, dim, OMEGA_H_CONSERVE)) return;
  auto kds2old_elems = old_mesh->ask_up(key_dim, dim);
  auto keys2old_elems = unmap_graph(keys2kds, kds2old_elems);
  auto keys2new_elems = Graph(keys2prods, prods2new_ents);
  std::map<Int, Graph> keys2old_elems_cat;
  std::map<Int, Graph> keys2new_elems_cat;
  if (old_mesh->has_tag(dim, "class_id")) {
    auto old_class_id = old_mesh->get_array<I32>(dim, "class_id");
    auto new_class_id = new_mesh->get_array<I32>(dim, "class_id");
    keys2old_elems_cat = categorize_graph(keys2old_elems, old_class_id);
    keys2new_elems_cat = categorize_graph(keys2new_elems, new_class_id);
  } else {
    keys2old_elems_cat[-1] = keys2old_elems;
    keys2new_elems_cat[-1] = keys2new_elems;
  }
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (tagbase->xfer() == OMEGA_H_CONSERVE) {
      transfer_conserve_tag(old_mesh, new_mesh, keys2old_elems_cat,
          keys2new_elems_cat, same_ents2old_ents, same_ents2new_ents, tagbase);
    }
  }
}

void fix_momentum_velocity_verts(Mesh* mesh, Int class_dim, I32 class_id,
    Int comp) {
  for (Int ent_dim = VERT; ent_dim <= mesh->dim(); ++ent_dim) {
    auto ent_marks = mark_class_closure(mesh, ent_dim, class_dim, class_id);
    auto comp_marks = multiply_each_by(I8(1 << comp), ent_marks);
    if (mesh->has_tag(ent_dim, "momentum_velocity_fixed")) {
      auto old_marks = mesh->get_array<I8>(ent_dim, "momentum_velocity_fixed");
      auto new_marks = bit_or_each(old_marks, comp_marks);
      mesh->set_tag(ent_dim, "momentum_velocity_fixed", new_marks);
    } else {
      mesh->add_tag(ent_dim, "momentum_velocity_fixed", 1, OMEGA_H_INHERIT,
          OMEGA_H_DO_OUTPUT, comp_marks);
    }
  }
}

bool has_fixed_momentum_velocity(Mesh* mesh) {
  return has_xfer(mesh, VERT, OMEGA_H_MOMENTUM_VELOCITY) &&
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
  auto edges2elems = mesh->ask_up(EDGE, mesh->dim());
  auto cands2elems = unmap_graph(cands2edges, edges2elems);
  auto cands2verts = get_closure_verts(mesh, cands2elems);
  auto vert_comps_are_fixed = mesh->get_array<I8>(VERT, "momentum_velocity_fixed");
  auto cand_vert_comps_are_fixed = unmap(cands2verts.ab2b, vert_comps_are_fixed, 1);
  auto cand_comps_are_fixed = fan_reduce_bit_and(cands2verts.a2ab,
      cand_vert_comps_are_fixed, 1);
  auto keep = each_eq_to(cand_comps_are_fixed, I8(0));
  return mesh->sync_subset_array(EDGE, keep, cands2edges, I8(0), 1);
}

template <Int dim>
class MomentumVelocity {
  Graph keys2target_verts;
  Graph keys2donor_verts;
  Graph target_verts2elems;
  Graph donor_verts2elems;
  Reals target_masses;
  Reals donor_masses;
  Reals donor_velocities;
  Write<Real> target_velocities;
  Read<I8> verts_are_fixed;

 public:
  MomentumVelocity(Mesh* donor_mesh, Mesh* target_mesh, Int key_dim,
      TagBase const* tagbase, LOs keys2kds, LOs keys2prods, LOs prods2new_elems,
      LOs same_verts2old_verts, LOs same_verts2new_verts) {
    CHECK(tagbase->ncomps() == dim);
    this->target_masses = target_mesh->get_array<Real>(dim, "mass");
    this->donor_masses = donor_mesh->get_array<Real>(dim, "mass");
    this->donor_velocities = to<Real>(tagbase)->array();
    this->target_velocities = Write<Real>(target_mesh->nverts() * dim);
    auto same_velocities = unmap(same_verts2old_verts, donor_velocities, dim);
    map_into(same_velocities, same_verts2new_verts, target_velocities, dim);
    auto keys2target_interior = Graph(keys2prods, prods2new_elems);
    this->keys2target_verts =
        get_closure_verts(target_mesh, keys2target_interior);
    auto kds2elems = donor_mesh->ask_up(key_dim, dim);
    auto keys2donor_interior = unmap_graph(keys2kds, kds2elems);
    this->keys2donor_verts = get_closure_verts(donor_mesh, keys2donor_interior);
    this->target_verts2elems = target_mesh->ask_up(VERT, dim);
    this->donor_verts2elems = donor_mesh->ask_up(VERT, dim);
    if (target_mesh->has_tag(VERT, "momentum_velocity_fixed")) {
      verts_are_fixed =
          target_mesh->get_array<I8>(VERT, "momentum_velocity_fixed");
    }
  }

  template <typename Arr>
  static DEVICE Vector<dim> get_interior_momentum(LO key,
      Graph const& keys2verts, Graph const& verts2elems, Reals const& masses,
      Arr const& velocities) {
    Vector<dim> momentum = zero_vector<dim>();
    for (auto kv = keys2verts.a2ab[key]; kv < keys2verts.a2ab[key + 1]; ++kv) {
      auto vert = keys2verts.ab2b[kv];
      auto velocity = get_vector<dim>(velocities, vert);
      for (auto ve = verts2elems.a2ab[vert]; ve < verts2elems.a2ab[vert + 1];
           ++ve) {
        auto elem = verts2elems.ab2b[ve];
        auto mass = masses[elem];
        momentum = momentum + mass * velocity / (dim + 1);
      }
    }
    return momentum;
  }

  DEVICE void operator()(LO key) const {
    auto donor_momentum = get_interior_momentum(key, keys2donor_verts,
        donor_verts2elems, donor_masses, donor_velocities);
    auto target_momentum = get_interior_momentum(key, keys2target_verts,
        target_verts2elems, target_masses, target_velocities);
    auto momentum_diff = (donor_momentum - target_momentum);
    auto begin = keys2target_verts.a2ab[key];
    auto end = keys2target_verts.a2ab[key + 1];
    Int ntarget_verts = 0;
    for (auto ktv = begin; ktv < end; ++ktv) {
      auto vert = keys2target_verts.ab2b[ktv];
      if (!verts_are_fixed.exists() || !verts_are_fixed[vert]) ++ntarget_verts;
    }
    CHECK(ntarget_verts > 0);
    auto momentum_factor = (dim + 1) * momentum_diff / ntarget_verts;
    for (auto ktv = begin; ktv < end; ++ktv) {
      auto vert = keys2target_verts.ab2b[ktv];
      if (verts_are_fixed.exists() && verts_are_fixed[vert]) continue;
      Real mass_sum = 0;
      for (auto ve = target_verts2elems.a2ab[vert];
           ve < target_verts2elems.a2ab[vert + 1]; ++ve) {
        auto elem = target_verts2elems.ab2b[ve];
        auto mass = target_masses[elem];
        mass_sum += mass;
      }
      auto velocity = get_vector<dim>(target_velocities, vert);
      velocity = velocity + (momentum_factor / mass_sum);
      set_vector(target_velocities, vert, velocity);
    }
  }

  static void apply(Mesh* donor_mesh, Mesh* target_mesh, Int key_dim,
      TagBase const* tagbase, LOs keys2kds, LOs keys2prods, LOs prods2new_elems,
      LOs same_verts2old_verts, LOs same_verts2new_verts) {
    MomentumVelocity<dim> self(donor_mesh, target_mesh, key_dim, tagbase,
        keys2kds, keys2prods, prods2new_elems, same_verts2old_verts,
        same_verts2new_verts);
    auto nkeys = keys2kds.size();
    parallel_for(nkeys, self);
    transfer_common3(target_mesh, VERT, tagbase, self.target_velocities);
  }
};

void transfer_momentum_velocity(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_elems, LOs same_verts2old_verts,
    LOs same_verts2new_verts) {
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tagbase = old_mesh->get_tag(VERT, i);
    if (tagbase->xfer() == OMEGA_H_MOMENTUM_VELOCITY) {
      if (dim == 3) {
        MomentumVelocity<3>::apply(old_mesh, new_mesh, key_dim, tagbase,
            keys2kds, keys2prods, prods2new_elems, same_verts2old_verts,
            same_verts2new_verts);
      }
      if (dim == 2) {
        MomentumVelocity<2>::apply(old_mesh, new_mesh, key_dim, tagbase,
            keys2kds, keys2prods, prods2new_elems, same_verts2old_verts,
            same_verts2new_verts);
      }
      CHECK(new_mesh->has_tag(VERT, tagbase->name()));
    }
  }
}

static Reals get_vertex_masses(Mesh* mesh) {
  CHECK(mesh->owners_have_all_upward(VERT));
  auto dim = mesh->dim();
  auto verts2elems = mesh->ask_up(VERT, dim);
  auto elems2mass = mesh->get_array<Real>(dim, "mass");
  auto verts2mass = graph_reduce(verts2elems, elems2mass, 1, OMEGA_H_SUM);
  verts2mass = multiply_each_by(1.0 / Real(dim), verts2mass);
  return mesh->sync_array(VERT, verts2mass, 1);
}

static Reals get_cavity_momenta(Mesh* mesh, Int key_dim, LOs keys2kds,
    Graph keys2elems, Reals vert_velocities) {
  auto dim = mesh->dim();
  auto elem_masses = mesh->get_array<Real>(dim, "mass");
  auto cavity_elem_masses = unmap(keys2elems.ab2b, elem_masses, 1);
  auto cavity_elem_velocities = average_field(mesh, dim, keys2elems.ab2b, dim,
      vert_velocities);
  auto cavity_elem_momenta = multiply_each(cavity_elem_velocities,
      cavity_elem_masses);
  auto cavity_momenta = fan_reduce(keys2elems.a2ab, cavity_elem_momenta, dim,
      OMEGA_H_SUM);
  return mesh->sync_subset_array(key_dim, cavity_momenta,
      keys2kds, 0.0, dim);
}

static void label_cavity_elems(Mesh* mesh, Graph keys2elems) {
  auto nlocal_keys = keys2elems.nnodes();
  auto offset = mesh->comm()->exscan(GO(nlocal_keys), OMEGA_H_SUM);
  auto key_globals = Read<GO>(nlocal_keys, offset, GO(1));
  auto elem_cavity_globals = map_onto(key_globals, keys2elems, mesh->nelems(),
      GO(-1), 1);
  mesh->add_tag(mesh->dim(), "cavity_id", 1, OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DO_OUTPUT, elem_cavity_globals);
}

static Read<I8> get_comps_are_fixed(Mesh* mesh) {
  if (mesh->has_tag(VERT, "momentum_velocity_fixed")) {
    return mesh->get_array<I8>(VERT, "momentum_velocity_fixed");
  } else {
    return Read<I8>(mesh->nverts(), I8(0));
  }
}

void do_momentum_velocity_elem_target(Mesh* donor_mesh, Mesh* target_mesh,
    Int key_dim, LOs keys2kds, LOs keys2prods, LOs prods2new_elems) {
  if (!has_xfer(donor_mesh, VERT, OMEGA_H_MOMENTUM_VELOCITY)) return;
  auto dim = donor_mesh->dim();
  auto nkeys = keys2kds.size();
  auto keys2target_elems = Graph(keys2prods, prods2new_elems);
  label_cavity_elems(target_mesh, keys2target_elems);
  auto kds2donor_elems = donor_mesh->ask_up(key_dim, dim);
  auto keys2donor_elems = unmap_graph(keys2kds, kds2donor_elems);
  auto keys2target_verts = get_closure_verts(target_mesh, keys2target_elems);
  auto comps_are_fixed = get_comps_are_fixed(target_mesh);
  for (Int tag_i = 0; tag_i < donor_mesh->ntags(VERT); ++tag_i) {
    auto tagbase = donor_mesh->get_tag(VERT, tag_i);
    if (tagbase->xfer() != OMEGA_H_MOMENTUM_VELOCITY) continue;
    CHECK(tagbase->ncomps() == dim);
    auto tag = to<Real>(tagbase);
    auto donor_vert_velocities = tag->array();
    auto donor_cavity_momenta = get_cavity_momenta(donor_mesh, key_dim,
        keys2kds, keys2donor_elems, donor_vert_velocities);
    auto target_vert_velocities = target_mesh->get_array<Real>(VERT, tag->name());
    auto target_cavity_momenta = get_cavity_momenta(donor_mesh, key_dim,
        keys2kds, keys2donor_elems, donor_vert_velocities);
    auto cavity_momentum_losses = subtract_each(donor_cavity_momenta,
        target_cavity_momenta);
    auto corrections_w = Write<Real>(target_mesh->nelems() * dim, 0.0);
    auto f = LAMBDA(LO key) {
      Few<Int, 3> nfree_verts;
      for (auto kv = keys2target_verts.a2ab[key];
           kv < keys2target_verts.a2ab[key + 1];
           ++kv) {
        auto v = keys2target_verts.ab2b[kv];
        auto code = Int(comps_are_fixed[v]);
        for (Int comp = 0; comp < dim; ++comp) {
          if (!((1 << comp) & code)) ++nfree_verts[comp];
        }
      }
      for (auto ke = keys2target_elems.a2ab[key];
           ke < keys2target_elems.a2ab[key + 1];
           ++ke) {
        auto e = keys2target_elems.ab2b[ke];
        for (Int comp = 0; comp < dim; ++comp) {
          corrections_w[e * dim + comp] =
            cavity_momentum_losses[key * dim + comp] / nfree_verts[comp];
        }
      }
    };
    parallel_for(nkeys, f);
    auto name = tag->name() + "_correction";
    target_mesh->add_tag(dim, name, dim, OMEGA_H_DONT_TRANSFER,
        OMEGA_H_DO_OUTPUT, Reals(corrections_w));
  }
}

void do_momentum_velocity_ghosted_target(Mesh* mesh) {
  if (!has_xfer(mesh, VERT, OMEGA_H_MOMENTUM_VELOCITY)) return;
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto dim = mesh->dim();
  auto cavity_ids = mesh->get_array<GO>(dim, "cavity_id");
  auto v2e = mesh->ask_up(VERT, dim);
  auto comps_are_fixed = get_comps_are_fixed(mesh);
  auto vert_masses = get_vertex_masses(mesh);
  for (Int tag_i = 0; tag_i < mesh->ntags(VERT); ++tag_i) {
    auto tagbase = mesh->get_tag(VERT, tag_i);
    if (tagbase->xfer() != OMEGA_H_MOMENTUM_VELOCITY) continue;
    CHECK(tagbase->ncomps() == dim);
    auto tag = to<Real>(tagbase);
    auto corr_name = tag->name() + "_correction";
    auto elem_corrections = mesh->get_array<Real>(dim, corr_name);
    auto marks = Write<I8>(mesh->nelems(), I8(0));
    auto vert_corrections_w = Write<Real>(mesh->nverts() * dim, 0.0);
    auto f = LAMBDA(LO v) {
      Vector<3> correction = vector_3(0.0, 0.0, 0.0);
      for (auto ve = v2e.a2ab[v]; ve < v2e.a2ab[v + 1]; ++ve) {
        auto e = v2e.ab2b[ve];
        if (marks[e]) continue;
        auto e_cav = cavity_ids[e];
        marks[e] = 1;
        for (auto ven = ve + 1; ven < v2e.a2ab[v + 1]; ++ven) {
          auto en = v2e.ab2b[ven];
          auto en_cav = cavity_ids[en];
          if (e_cav == en_cav) marks[en] = 1;
        }
        for (Int comp = 0; comp < dim; ++comp) {
          correction[comp] += elem_corrections[e * dim + comp];
        }
      }
      for (Int comp = 0; comp < dim; ++comp) {
        if (!comps_are_fixed[v * dim + comp]) {
          vert_corrections_w[v * dim + comp] = correction[comp];
        }
      }
    };
    parallel_for(mesh->nverts(), f);
    auto vert_corrections = Reals(vert_corrections_w);
    vert_corrections = mesh->sync_array(dim, vert_corrections, dim);
    auto old_velocities = tag->array();
    auto velocity_corrections = divide_each(vert_corrections, vert_masses);
    auto new_velocities = add_each(old_velocities, velocity_corrections);
    mesh->set_tag(VERT, tag->name(), new_velocities);
    mesh->remove_tag(dim, corr_name);
  }
  mesh->remove_tag(dim, "cavity_id");
}

}  // end namespace Omega_h
