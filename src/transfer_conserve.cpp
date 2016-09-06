#include "transfer_conserve.hpp"

#include "Omega_h_r3d.hpp"
#include "loop.hpp"
#include "size.hpp"
#include "tag.hpp"
#include "transfer.hpp"
#include "indset.hpp"

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
static void transfer_conserve_tmpl(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents, TagBase const* tagbase) {
  auto name = tagbase->name();
  auto old_tag = to<Real>(tagbase);
  auto ncomps = old_tag->ncomps();
  auto old_data = old_tag->array();
  auto kds2elems = old_mesh->ask_up(key_dim, dim);
  auto kds2kd_elems = kds2elems.a2ab;
  auto kd_elems2elems = kds2elems.ab2b;
  auto measure = RealElementSizes(new_mesh);
  auto elem_verts2verts = new_mesh->ask_verts_of(dim);
  auto nkeys = keys2kds.size();
  auto nprods = keys2prods.last();
  auto prod_data_w = Write<Real>(nprods * ncomps);
  auto f = LAMBDA(LO key) {
    auto kd = keys2kds[key];
    Real total_new_size = 0.0;
    for (auto prod = keys2prods[key]; prod < keys2prods[key + 1]; ++prod) {
      auto new_elem = prods2new_ents[prod];
      auto v = gather_verts<dim + 1>(elem_verts2verts, new_elem);
      auto size = measure.measure(v);
      total_new_size += size;
    }
    for (Int comp = 0; comp < ncomps; ++comp) {
      Real sum = 0.0;
      for (auto kd_elem = kds2kd_elems[kd]; kd_elem < kds2kd_elems[kd + 1];
           ++kd_elem) {
        auto old_elem = kd_elems2elems[kd_elem];
        sum += old_data[old_elem * ncomps + comp];
      }
      for (auto prod = keys2prods[key]; prod < keys2prods[key + 1]; ++prod) {
        auto new_elem = prods2new_ents[prod];
        auto v = gather_verts<dim + 1>(elem_verts2verts, new_elem);
        auto size = measure.measure(v);
        prod_data_w[prod * ncomps + comp] = sum * (size / total_new_size);
      }
    }
  };
  parallel_for(nkeys, f);
  auto prod_data = Reals(prod_data_w);
  transfer_common(old_mesh, new_mesh, dim, same_ents2old_ents,
      same_ents2new_ents, prods2new_ents, old_tag, prod_data);
}

void transfer_conserve(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  auto dim = new_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (tagbase->xfer() == OMEGA_H_CONSERVE) {
      if (dim == 3) {
        transfer_conserve_tmpl<3>(old_mesh, new_mesh, key_dim, keys2kds,
            keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
            tagbase);
      } else if (dim == 2) {
        transfer_conserve_tmpl<2>(old_mesh, new_mesh, key_dim, keys2kds,
            keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
            tagbase);
      }
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

Read<I8> find_buffered_indset(
    Mesh* mesh, Int key_dim,
    Reals qualities,
    Read<I8> unbuffered_indset) {
  auto unbuffered_conflicts = mesh->ask_star(key_dim);
  auto buffered_conflicts = get_buffered_conflict_graph(unbuffered_conflicts,
      unbuffered_indset);
  return find_indset(mesh, key_dim, buffered_conflicts,
      qualities, unbuffered_indset);
}

}  // end namespace Omega_h
