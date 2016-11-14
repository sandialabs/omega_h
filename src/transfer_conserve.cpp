#include "transfer_conserve.hpp"

#include "Omega_h_math.hpp"
#include "Omega_h_r3d.hpp"
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

bool needs_buffer_layers(Mesh* mesh) {
  for (Int i = 0; i < mesh->ntags(VERT); ++i)
    if (mesh->get_tag(VERT, i)->xfer() == OMEGA_H_MOMENTUM_VELOCITY)
      return true;
  return false;
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
  auto dim = mesh->dim();
  auto all_comps = I8((1 << dim) - 1);
  auto out = Write<I8>(ncands);
  auto f = LAMBDA(LO cand) {
    auto e = cands2edges[cand];
    auto code = cand_codes[cand];
    for (Int eev_col = 0; eev_col < 2; ++eev_col) {
      if (!collapses(code, eev_col)) continue;
      auto v_col = ev2v[e * 2 + eev_col];
      I8 free_comps = 0;
      for (auto ve = v2e.a2ab[v_col]; ve < v2e.a2ab[v_col + 1]; ++ve) {
        auto e2 = v2e.ab2b[ve];
        auto e2_code = v2e.codes[ve];
        auto eev_in = code_which_down(e2_code);
        auto eev_out = 1 - eev_in;
        auto ov = ev2v[e2 * 2 + eev_out];
        free_comps |= ~(comps_are_fixed[ov]);
      }
      if (free_comps & all_comps != all_comps) {
        code = dont_collapse(code, eev_col);
      }
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
  auto vert_comps_are_free = bit_neg_each(vert_comps_are_fixed);
  auto cand_comps_are_free = graph_reduce(cands2verts, vert_comps_are_fixed, 1,
      OMEGA_H_BIT_OR);
  auto keep = each_eq_to(cand_comps_are_free, ~I8(0));
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
  auto dim = donor_mesh->dim();
  auto verts2elems = donor_mesh->ask_up(VERT, dim);
  auto elems2mass = donor_mesh->get_array<Real>(dim, "mass");
  auto verts2mass = graph_reduce(verts2elems, elems2mass, 1, OMEGA_H_SUM);
  verts2mass = multiply_each_by(1.0 / Real(dim), verts2mass);
  return mesh->sync_array(VERT, verts2mass, 1);
}

void do_momentum_velocity_ghosted_donor(Mesh* donor_mesh) {
  if (!has_xfer(donor_mesh, VERT, OMEGA_H_MOMENTUM_VELOCITY)) return;
  auto dim = donor_mesh->dim();
  auto vert_masses = get_vertex_masses(donor_mesh);
  for (Int i = 0; i < donor_mesh->ntags(VERT); ++i) {
    auto tagbase = donor_mesh->get_tag(VERT, i);
    if (tagbase->xfer() != OMEGA_H_MOMENTUM_VELOCITY) continue;
    auto tag = to<Real>(tagbase);
    auto vert_velocities = tag->array();
    auto vert_momenta = multiply_each(vert_velocities, vert_masses);
    auto name = tagbase->name() + "_momentum";
    donor_mesh->add_tag(VERT, name, dim, OMEGA_H_LINEAR_INTERP,
        tagbase->outflag(), vert_momenta);
  }
}

static Reals get_cavity_momenta(Mesh* mesh, Graph keys2elems,
    Reals vert_velocities) {
  auto dim = mesh->dim();
  auto elem_masses = mesh->get_array<Real>(dim, "mass");
  auto cavity_elem_velocities = average_field(mesh, dim, keys2elems.ab2b, dim,
      vert_velocities);
  auto cavity_elem_momenta = multiply_each(elem_velocities, elem_masses);
  return fan_reduce(keys2elems.a2ab, cavity_elem_momenta, dim,
      OMEGA_H_SUM);
}

void do_momentum_velocity_elem_target(Mesh* donor_mesh, Mesh* target_mesh,
    Int key_dim, LOs keys2kds, LOs keys2prods, LOs prods2new_elems,
    LOs same_verts2old_verts, LOs same_verts2new_verts) {
  if (!has_xfer(donor_mesh, VERT, OMEGA_H_MOMENTUM_VELOCITY)) return;
  auto dim = donor_mesh->dim();
  auto donor_masses = donor_mesh->get_array<Real>(dim, "mass");
  auto target_masses = target_mesh->get_array<Real>(dim, "mass");
  auto keys2target_elems = Graph(keys2prods, prods2new_elems);
  auto keys2target_verts = get_closure_verts(target_mesh, keys2target_elems);
  auto target_elems2verts = target_mesh->ask_verts_of_elems();
  auto kds2elems = donor_mesh->ask_up(key_dim, dim);
  auto keys2donor_elems = unmap_graph(keys2kds, kds2elems);
  auto keys2donor_verts = get_closure_verts(donor_mesh, keys2donor_elems);
  auto donor_elems2verts = donor_mesh->ask_verts_of_elems();
  for (Int tag_i = 0; tag_i < donor_mesh->ntags(VERT); ++tag_i) {
    auto tagbase = donor_mesh->get_tag(VERT, tag_i);
    if (tagbase->xfer() != OMEGA_H_MOMENTUM_VELOCITY) continue;
    auto tag = to<Real>(tagbase);
    auto donor_velocities = tag->array();
    auto target_velocities = target_mesh->get_array<Real>(VERT, tag->name());
  }
}

}  // end namespace Omega_h
