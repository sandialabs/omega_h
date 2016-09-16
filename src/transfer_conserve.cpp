#include "transfer_conserve.hpp"

#include "Omega_h_r3d.hpp"
#include "loop.hpp"
#include "size.hpp"
#include "tag.hpp"
#include "transfer.hpp"
#include "map.hpp"
#include "graph.hpp"
#include "simplices.hpp"

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
static void transfer_conserve_dim(
    Mesh* new_mesh,
    TagBase const* tagbase,
    Graph keys2old_mat_elems,
    Graph keys2new_mat_elems,
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
         kne < keys2new_mat_elems.a2ab[key + 1];
         ++kne) {
      auto new_elem = keys2new_mat_elems.ab2b[kne];
      auto v = gather_verts<dim + 1>(elem_verts2verts, new_elem);
      auto size = measure.measure(v);
      total_new_size += size;
    }
    for (Int comp = 0; comp < ncomps; ++comp) {
      Real sum = 0.0;
      for (auto koe = keys2old_mat_elems.a2ab[key];
           koe < keys2old_mat_elems.a2ab[key + 1];
           ++koe) {
        auto old_elem = keys2old_mat_elems.ab2b[koe];
        sum += old_data[old_elem * ncomps + comp];
      }
      for (auto kne = keys2new_mat_elems.a2ab[key];
           kne < keys2new_mat_elems.a2ab[key + 1];
           ++kne) {
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
    std::map<Int, Graph> keys2new_elems_cat,
    LOs same_ents2old_ents, LOs same_ents2new_ents,
    TagBase const* tagbase) {
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
      transfer_conserve_dim<3>(new_mesh, tagbase,
          keys2old_mat_elems, keys2new_mat_elems, new_data_w);
    if (dim == 2)
      transfer_conserve_dim<2>(new_mesh, tagbase,
          keys2old_mat_elems, keys2new_mat_elems, new_data_w);
  }
  transfer_common2(old_mesh, new_mesh, dim, same_ents2old_ents,
      same_ents2new_ents, tagbase, new_data_w);
}

static bool has_xfer(Mesh* mesh, Int dim, Omega_h_Xfer xfer) {
  for (Int i = 0; i < mesh->ntags(dim); ++i)
    if (mesh->get_tag(dim, i)->xfer() == xfer)
      return true;
  return false;
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
      transfer_conserve_tag(old_mesh, new_mesh,
          keys2old_elems_cat, keys2new_elems_cat,
          same_ents2old_ents, same_ents2new_ents, tagbase);
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

static bool starts_with(std::string const& a, std::string const& b) {
  return 0 == a.compare(0, b.length(), b);
}

static std::string remove_prefix(std::string const& a, std::string const& prefix) {
  return a.substr(prefix.length(), std::string::npos);
}

class MomentumVelocity {
protected:
  LOs target_elem_verts2verts;
  LOs donor_elem_verts2verts;
  LOs verts2dofs;
  Graph keys2target_buffer;
  Graph keys2target_interior;
  Graph keys2donor_interior;
  LOs target_elems2donor_elems;
  Reals target_coords;
  Reals donor_coords;

public:
  MomentumVelocity(Mesh* donor_mesh, Mesh* target_mesh,
      Int key_dim, LOs keys2kds,
      LOs keys2prods, LOs prods2new_ents,
      LOs same_ents2old_ents, LOs same_ents2new_ents) {
    auto elem_dim = donor_mesh->dim();
    this->target_elem_verts2verts = target_mesh->ask_verts_of(elem_dim);
    this->donor_elem_verts2verts = donor_mesh->ask_verts_of(elem_dim);
    this->keys2target_interior = Graph(keys2prods, prods2new_ents);
    auto keys2target_verts = get_closure_verts(target_mesh, keys2target_interior);
    this->verts2dofs = number_cavity_ents(target_mesh, keys2target_verts,
        elem_dim);
    auto nkeys = keys2kds.size();
    auto nkds = donor_mesh->nents(key_dim);
    auto kds_are_keys = map_onto<I8>(Read<I8>(nkeys, 1), keys2kds, nkds, 0, 1);
    auto kds2donor_elems = get_buffered_elems(donor_mesh, key_dim, kds_are_keys);
    auto keys2donor_elems = unmap_graph(keys2kds, kds2donor_elems);
    auto ndonor_elems = donor_mesh->nelems();
    auto donor_elems2target_elems = map_onto<LO>(same_ents2new_ents,
      same_ents2old_ents, ndonor_elems, -1, 1);
    this->keys2target_buffer = get_target_buffer_elems(keys2donor_elems,
        donor_elems2target_elems);
    auto ntarget_elems = target_mesh->nelems();
    this->target_elems2donor_elems = map_onto<LO>(same_ents2old_ents,
      same_ents2new_ents, ntarget_elems, -1, 1);
    this->target_coords = target_mesh->coords();
    this->donor_coords = donor_mesh->coords();
  };
};

template <Int dim>
class MomentumVelocityDim : public MomentumVelocity {
public:
  static constexpr Int nverts_per_elem = dim + 1;
  static constexpr Int max_dofs = (AvgDegree<dim,0,1>::value + 1) * 2;
  static constexpr Real coupling_factor =
    ParentElementSize<dim>::value / (dim + 1) * (dim + 2);
  using MassMatrix = Matrix<max_dofs, max_dofs>;
  using RHS = Matrix<max_dofs, dim>;

protected:
  Reals donor_velocities;
  Reals donor_densities;
  Reals target_densities;
  Reals target_sizes;

public:
  MomentumVelocityDim(
      MomentumVelocity parent,
      Mesh* donor_mesh,
      Mesh* target_mesh,
      TagBase const* tagbase):
      MomentumVelocity(parent) {
    auto velocity_name = tagbase->name();
    if (!starts_with(velocity_name, "velocity")) {
      Omega_h_fail("%s tranferred as momentum-velocity,"
                   " but name needs to start with \"velocity\"\n",
                   velocity_name.c_str());
    }
    auto suffix = remove_prefix(velocity_name, "velocity");
    auto density_name = std::string("density") + suffix;
    this->donor_velocities = donor_mesh->get_array<Real>(VERT, velocity_name);
    this->donor_densities = donor_mesh->get_array<Real>(dim, density_name);
    this->target_densities = target_mesh->get_array<Real>(dim, density_name);
    this->target_sizes = measure_elements_real(target_mesh);
  }

protected:

  DEVICE void elem_into_mass_matrix(LO elem, MassMatrix& A) {
    auto rho_a = target_densities[elem];
    auto V_a = target_sizes[elem];
    auto contrib = coupling_factor * rho_a * V_a;
    for (Int elem_vert1 = 0; elem_vert1 < nverts_per_elem; ++elem_vert1) {
      auto vert1 = target_elem_verts2verts[elem * nverts_per_elem + elem_vert1];
      auto dof1 = verts2dofs[vert1];
      if (dof1 < 0) continue;
      for (Int elem_vert2 = 0; elem_vert2 < nverts_per_elem;
           ++elem_vert2) {
        if (elem_vert2 == elem_vert1) continue;
        auto vert2 = target_elem_verts2verts[elem * nverts_per_elem + elem_vert2];
        auto dof2 = verts2dofs[vert2];
        if (dof2 < 0) continue;
        A[dof1][dof2] += contrib;
      }
    }
  }

  DEVICE void buffer_elem_into_rhs(LO target_elem, RHS& b) {
    auto rho_a = target_densities[target_elem];
    auto V_a = target_sizes[target_elem];
    auto contrib = coupling_factor * rho_a * V_a;
    auto donor_elem = target_elems2donor_elems[target_elem];
    for (Int elem_vert1 = 0; elem_vert1 < nverts_per_elem; ++elem_vert1) {
      auto target_vert1 = target_elem_verts2verts[
        target_elem * nverts_per_elem + elem_vert1];
      auto dof1 = verts2dofs[target_vert1];
      if (dof1 < 0) continue;
      auto donor_vert = donor_elem_verts2verts[
        donor_elem * nverts_per_elem + elem_vert1];
      auto donor_velocity = get_vector<dim>(donor_velocities, donor_vert);
      for (Int elem_vert2 = 0; elem_vert2 < nverts_per_elem;
           ++elem_vert2) {
        if (elem_vert2 == elem_vert1) continue;
        auto target_vert2 = target_elem_verts2verts[
          target_elem * nverts_per_elem + elem_vert2];
        auto dof2 = verts2dofs[target_vert2];
        if (dof2 < 0) continue;
        for (Int i = 0; i < dim; ++i)
          b[i][dof2] += donor_velocity[i] * contrib;
      }
    }
  }

  DEVICE void elem_pair_into_rhs(LO target_elem, LO donor_elem, RHS& b) {
    auto target_pts = gather_vectors<dim + 1, dim>(
        target_coords, target_elem_verts2verts);
    auto donor_pts = gather_vectors<dim + 1, dim>(
        donor_coords, donor_elem_verts2verts);
    auto domain = r3d::intersect_simplices(target_pts, donor_pts);
    for (Int elem_vert1 = 0; elem_vert1 < nverts_per_elem; ++elem_vert1) {
      for (Int elem_vert2 = 0; elem_vert2 < nverts_per_elem;
           ++elem_vert2) {
        if (elem_vert2 == elem_vert1) continue;
        /* TODO ! */
      }
    }
  }

  DEVICE MassMatrix elems_into_mass_matrix(LO key, Graph const& keys2elems,
      MassMatrix& A) {
    for (auto ke = keys2elems.a2ab[key];
         ke < keys2elems.a2ab[key + 1];
         ++ke) {
      auto target_elem = keys2elems.ab2b[ke];
      elem_into_mass_matrix(target_elem, A);
    }
    return A;
  }

  DEVICE MassMatrix assemble_mass_matrix(LO key) {
    MassMatrix A = zero_matrix<max_dofs, max_dofs>();
    elems_into_mass_matrix(key, keys2target_interior, A);
    elems_into_mass_matrix(key, keys2target_buffer, A);
    return A;
  }

};

}  // end namespace Omega_h
