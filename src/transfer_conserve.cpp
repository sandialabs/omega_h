#include "transfer_conserve.hpp"

#include "Omega_h_r3d.hpp"
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
static void transfer_conserve_dim(Mesh* new_mesh, TagBase const* tagbase,
    Graph keys2old_mat_elems, Graph keys2new_mat_elems,
    Write<Real> new_data_w) {
  auto ncomps = tagbase->ncomps();
  auto old_tag = to<Real>(tagbase);
  auto old_data = old_tag->array();
  auto elem_verts2verts = new_mesh->ask_verts_of(dim);
  auto measure = RealElementSizes(new_mesh);
  auto nkeys = keys2old_mat_elems.nnodes();
  auto f = LAMBDA(LO key) {
    Real total_new_size = 0.0;
    for (auto kne = keys2new_mat_elems.a2ab[key];
         kne < keys2new_mat_elems.a2ab[key + 1]; ++kne) {
      auto new_elem = keys2new_mat_elems.ab2b[kne];
      auto v = gather_verts<dim + 1>(elem_verts2verts, new_elem);
      auto size = measure.measure(v);
      total_new_size += size;
    }
    for (Int comp = 0; comp < ncomps; ++comp) {
      Real sum = 0.0;
      for (auto koe = keys2old_mat_elems.a2ab[key];
           koe < keys2old_mat_elems.a2ab[key + 1]; ++koe) {
        auto old_elem = keys2old_mat_elems.ab2b[koe];
        sum += old_data[old_elem * ncomps + comp];
      }
      for (auto kne = keys2new_mat_elems.a2ab[key];
           kne < keys2new_mat_elems.a2ab[key + 1]; ++kne) {
        auto new_elem = keys2new_mat_elems.ab2b[kne];
        auto v = gather_verts<dim + 1>(elem_verts2verts, new_elem);
        auto size = measure.measure(v);
        new_data_w[new_elem * ncomps + comp] = sum * (size / total_new_size);
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
      transfer_conserve_dim<3>(new_mesh, tagbase, keys2old_mat_elems,
          keys2new_mat_elems, new_data_w);
    if (dim == 2)
      transfer_conserve_dim<2>(new_mesh, tagbase, keys2old_mat_elems,
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

template <Int dim>
static void transfer_conserve_r3d_tmpl(Mesh* old_mesh, Mesh* new_mesh,
    Int key_dim, LOs keys2kds, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents, TagBase const* tagbase) {
  CHECK(old_mesh->dim() == dim);
  auto name = tagbase->name();
  auto old_tag = to<Real>(tagbase);
  auto ncomps = old_tag->ncomps();
  auto old_data = old_tag->array();
  auto kds2elems = old_mesh->ask_up(key_dim, dim);
  auto kds2kd_elems = kds2elems.a2ab;
  auto kd_elems2elems = kds2elems.ab2b;
  auto old_ev2v = old_mesh->ask_elem_verts();
  auto old_coords = old_mesh->coords();
  auto new_ev2v = new_mesh->ask_elem_verts();
  auto new_coords = new_mesh->coords();
  auto nkeys = keys2kds.size();
  auto nprods = keys2prods.last();
  auto prod_data_w = Write<Real>(nprods * ncomps);
  auto f = LAMBDA(LO key) {
    auto kd = keys2kds[key];
    for (auto prod = keys2prods[key]; prod < keys2prods[key + 1]; ++prod) {
      auto target_elem = prods2new_ents[prod];
      auto target_verts = gather_verts<dim + 1>(new_ev2v, target_elem);
      auto target_points =
          gather_vectors<dim + 1, dim>(new_coords, target_verts);
      auto target_p = &prod_data_w[prod * ncomps];
      for (Int comp = 0; comp < ncomps; ++comp) target_p[comp] = 0;
      /* total_size will track the sum of the volumes of donor element
       * intersections with this target element.
       * if the cavity boundaries are the same, it should end up equal
       * to target element volume, but if this is an edge collapse on
       * a curved boundary it could be less and we use it as the denominator
       * to avoid losing conservation in that case.
       */
      Real total_size = 0;
      for (auto kd_elem = kds2kd_elems[kd]; kd_elem < kds2kd_elems[kd + 1];
           ++kd_elem) {
        auto donor_elem = kd_elems2elems[kd_elem];
        auto donor_p = &old_data[donor_elem * ncomps];
        auto donor_verts = gather_verts<dim + 1>(old_ev2v, donor_elem);
        auto donor_points =
            gather_vectors<dim + 1, dim>(old_coords, donor_verts);
        auto intersection =
            r3d::intersect_simplices(target_points, donor_points);
        auto intersection_size = r3d::measure(intersection);
        for (Int comp = 0; comp < ncomps; ++comp) {
          target_p[comp] += intersection_size * donor_p[comp];
        }
        total_size += intersection_size;
      }
      for (Int comp = 0; comp < ncomps; ++comp) target_p[comp] /= total_size;
    }
  };
  parallel_for(nkeys, f);
  auto prod_data = Reals(prod_data_w);
  transfer_common(old_mesh, new_mesh, dim, same_ents2old_ents,
      same_ents2new_ents, prods2new_ents, old_tag, prod_data);
}

void transfer_conserve_r3d(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  auto dim = new_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (tagbase->xfer() == OMEGA_H_CONSERVE_R3D) {
      if (dim == 3) {
        transfer_conserve_r3d_tmpl<3>(old_mesh, new_mesh, key_dim, keys2kds,
            keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
            tagbase);
      } else if (dim == 2) {
        transfer_conserve_r3d_tmpl<2>(old_mesh, new_mesh, key_dim, keys2kds,
            keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
            tagbase);
      }
    }
  }
}

void transfer_conserve_r3d_refine(Mesh* old_mesh, Mesh* new_mesh,
    LOs keys2edges, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (tagbase->xfer() == OMEGA_H_CONSERVE_R3D) {
      transfer_inherit_refine<Real>(old_mesh, new_mesh, keys2edges, dim,
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase->name());
    }
  }
}

bool needs_buffer_layers(Mesh* mesh) {
  for (Int i = 0; i < mesh->ntags(VERT); ++i)
    if (mesh->get_tag(VERT, i)->xfer() == OMEGA_H_MOMENTUM_VELOCITY)
      return true;
  return false;
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

 public:
  MomentumVelocity(Mesh* donor_mesh, Mesh* target_mesh,
      Int key_dim, TagBase const* tagbase,
      LOs keys2kds, LOs keys2prods, LOs prods2new_elems,
      LOs same_verts2old_verts, LOs same_verts2new_verts) {
    auto velocity_name = tagbase->name();
    if (velocity_name != "velocity") {
      Omega_h_fail(
          "%s tranferred as momentum-conserving velocity,"
          " but its name is not \"velocity\"\n",
          velocity_name);
    }
    this->target_masses = target_mesh->get_array<Real>(dim, "mass");
    this->donor_masses = donor_mesh->get_array<Real>(dim, "mass");
    this->donor_velocities = donor_mesh->get_array<Real>(VERT, velocity_name);
    this->target_velocities = Write<Real>(target_mesh->nverts() * dim);
    auto same_velocities = unmap(same_verts2old_verts, donor_velocities, dim);
    map_into(same_velocities, same_verts2new_verts, target_velocities, dim);
    auto keys2target_interior = Graph(keys2prods, prods2new_elems);
    this->keys2target_verts =
        get_closure_verts(target_mesh, keys2target_interior);
    auto kds2elems = donor_mesh->ask_up(key_dim, dim);
    auto keys2donor_interior = unmap_graph(keys2kds, kds2elems);
    this->keys2donor_verts =
        get_closure_verts(donor_mesh, keys2donor_interior);
    this->target_verts2elems = target_mesh->ask_up(VERT, dim);
    this->donor_verts2elems = donor_mesh->ask_up(VERT, dim);
  }

  template <typename Arr>
  static DEVICE Vector<dim> get_interior_momentum(
      LO key,
      Graph const& keys2verts,
      Graph const& verts2elems,
      Reals const& masses,
      Arr const& velocities) {
    Vector<dim> momentum = zero_vector<dim>();
    for (auto kv = keys2verts.a2ab[key];
         kv < keys2verts.a2ab[key + 1];
         ++kv) {
      auto vert = keys2verts.ab2b[kv];
      auto velocity = get_vector<dim>(velocities, vert);
      for (auto ve = verts2elems.a2ab[vert];
           ve < verts2elems.a2ab[vert + 1];
           ++ve) {
        auto elem = verts2elems.ab2b[ve];
        auto mass = masses[elem];
        momentum = momentum + mass * velocity / 4.0;
      }
    }
    return momentum;
  }

  DEVICE void operator()(LO key) const {
    auto donor_momentum = get_interior_momentum(
        key, keys2donor_verts, donor_verts2elems,
        donor_masses, donor_velocities);
    auto target_momentum = get_interior_momentum(
        key, keys2target_verts, target_verts2elems,
        target_masses, target_velocities);
    for (auto kv = keys2target_verts.a2ab[key];
         kv < keys2target_verts.a2ab[key + 1];
         ++kv) {
      auto vert = keys2target_verts.ab2b[kv];
      auto velocity = get_vector<dim>(target_velocities, vert);
      for (Int i = 0; i < dim; ++i) {
        velocity[i] = velocity[i] * target_momentum[i] / donor_momentum[i];
      }
      set_vector<dim>(target_velocities, vert, velocity);
    }
  }

  static void apply(Mesh* donor_mesh, Mesh* target_mesh,
      Int key_dim, TagBase const* tagbase,
      LOs keys2kds, LOs keys2prods, LOs prods2new_elems,
      LOs same_verts2old_verts, LOs same_verts2new_verts) {
    MomentumVelocity<dim> self(donor_mesh, target_mesh, key_dim, tagbase,
        keys2kds, keys2prods, prods2new_elems,
        same_verts2old_verts, same_verts2new_verts);
    auto nkeys = keys2kds.size();
    parallel_for(nkeys, self);
    transfer_common3(target_mesh, VERT, tagbase, self.target_velocities);
  }

};

void transfer_momentum_velocity(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_elems,
    LOs same_verts2old_verts, LOs same_verts2new_verts) {
  auto dim = new_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (tagbase->xfer() == OMEGA_H_MOMENTUM_VELOCITY) {
      if (dim == 3) {
        MomentumVelocity<3>::apply(old_mesh, new_mesh, key_dim, tagbase,
            keys2kds, keys2prods, prods2new_elems,
            same_verts2old_verts, same_verts2new_verts);
      }
      if (dim == 2) {
        MomentumVelocity<2>::apply(old_mesh, new_mesh, key_dim, tagbase,
            keys2kds, keys2prods, prods2new_elems,
            same_verts2old_verts, same_verts2new_verts);
      }
    }
  }
}

}  // end namespace Omega_h
