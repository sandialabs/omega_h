#include "Omega_h_conserve.hpp"

#include "Omega_h_adj.hpp"
#include "Omega_h_graph.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_r3d.hpp"

namespace Omega_h {

/* computes the error in the integral over a cavity, given the old
   and new densities at elements.
   also factors in any previous recorded error on old elements */
static void track_subcavs_integral_error(Mesh* old_mesh, Mesh* new_mesh,
    Graph keys2old_elems, Graph keys2new_elems,
    Reals old_elem_densities, Reals new_elem_densities, std::string const& error_name,
    Write<Real>* new_elem_errors_w, bool conserves_integrals) {
  Reals cav_errors;
  if (!conserves_integrals) {
    auto dim = old_mesh->dim();
    auto ncomps = old_elem_densities.size() / old_mesh->nelems();
    auto old_elem_sizes = old_mesh->ask_sizes();
    auto new_elem_sizes = new_mesh->ask_sizes();
    auto old_cav_elem_densities = unmap(
        keys2old_elems.ab2b, old_elem_densities, ncomps);
    auto new_cav_elem_densities = unmap(
        keys2new_elems.ab2b, old_elem_densities, ncomps);
    auto old_cav_elem_sizes = unmap(keys2old_elems.ab2b, old_elem_sizes, 1);
    auto new_cav_elem_sizes = unmap(keys2new_elems.ab2b, old_elem_sizes, 1);
    auto old_cav_elem_integrals = multiply_each(
        old_cav_elem_densities, old_cav_elem_sizes);
    auto new_cav_elem_integrals = multiply_each(
        new_cav_elem_densities, new_cav_elem_sizes);
    auto old_cav_integrals = fan_reduce(
        keys2old_elems.a2ab, old_cav_elem_integrals, ncomps, OMEGA_H_SUM);
    auto new_cav_integrals = fan_reduce(
        keys2new_elems.a2ab, new_cav_elem_integrals, ncomps, OMEGA_H_SUM);
    cav_errors = subtract_each(new_cav_integrals, old_cav_integrals);
  }
  if (old_mesh->has_tag(dim, error_name)) {
    auto old_elem_errors = old_mesh->get_array<Real>(dim, error_name);
    auto old_cav_errors = graph_reduce(
        keys2old_elems, old_elem_errors, ncomps, OMEGA_H_SUM);
    if (!conserves_integrals) cav_errors = add_each(cav_errors, old_cav_errors);
    else cav_errors = old_cav_errors;
  }
  if (!cav_errors.exists()) return;
  auto new_cav_sizes = fan_reduce(
      keys2new_elems.a2ab, new_cav_elem_sizes, 1, OMEGA_H_SUM);
  auto new_cav_error_densities = divide_each(cav_errors, new_cav_sizes);
  auto new_cav_elem_error_densities = expand(new_cav_error_densities,
      keys2new_elems, ncomps);
  auto new_cav_elem_errors = multiply_each(
      new_cav_elem_error_densities, new_cav_elem_sizes);
  if (!new_elem_errors_w.exists()) {
    *new_elem_errors_w = Write<Real>(new_mesh->nelems() * ncomps, 0.0);
  }
  map_into(new_cav_elem_errors, keys2new_elems.ab2b, *new_elem_errors_w, ncomps);
}

static void track_cavs_integral_error(Mesh* old_mesh, Mesh* new_mesh,
    std::vector<std::pair<Graph, Graph>> mats_keys2elems,
    Reals old_elem_densities, Reals new_elem_densities, std::string const& error_name,
    LOS same_ents2old_ents, LOS same_ents2new_ents,
    bool conserves_integrals) {
  Write<Real> new_elem_errors_w;
  if (old_mesh->has_tag(dim, error_name)) {
    auto old_tag = mesh->get_tag<Real>(dim, error_name);
    auto ncomps = old_tag->ncomps;
    new_elem_errors_w = Write<Real>(new_mesh->nelems() * ncomps);
    auto old_elem_errors = old_tag->array();
    auto same_errors = unmap(same_ents2old_ents, old_elem_errors);
    map_into(same_errors, same_ents2new_ents, new_elem_errors_w, ncomps);
  }
  for (auto pair : mats_keys2elems) {
    track_subcavs_integral_error(old_mesh, new_mesh, pair.first, pair.second,
        old_elem_densities, new_elem_densities, error_name, &new_elem_errors_w,
        conserves_integrals);
  }
  if (new_elem_errors_w.exists()) {
    new_mesh->add_tag(dim, error_name, 1, new_elem_errors_w);
  }
}

static void track_cavs_density_error(
    Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
    std::vector<std::pair<Graph, Graph>> mats_keys2elems, TagBase const* tagbase,
    LOs same_ents2old_ents, LOs same_ents2new_ents, bool conserves_density) {
  auto dim = old_mesh->dim();
  auto ncomps = tagbase->ncomps();
  auto integral_name = opts.integral_map.find(tagbase->name()).second;
  auto error_name = integral_name + "_error";
  auto old_tag = to<Real>(tagbase);
  auto old_elem_densities = old_tag->array();
  auto new_elem_densities = new_mesh->get_array<Real>(dim, tagbase->name());
  track_cavs_integral_error(old_mesh, new_mesh, mats_keys2elems,
      old_elem_densities, new_elem_densities, error_name, conserves_density);
}

static void track_cavs_size_error(Mesh* old_mesh, Mesh* new_mesh,
    std::vector<std::pair<Graph, Graph>> mats_keys2elems,
    LOs same_ents2old_ents, LOs same_ents2new_ents, bool conserves_size) {
  auto dim = old_mesh->dim();
  auto ncomps = tagbase->ncomps();
  auto error_name = "size_error";
  auto old_elem_densities = Reals(old_mesh->nelems(), 1.0);
  auto new_elem_densities = Reals(new_mesh->nelems(), 1.0);
  track_cavs_integral_error(old_mesh, new_mesh, mats_keys2elems,
      old_elem_densities, new_elem_densities, error_name, conserves_size);
}

static void track_cavs_momentum_error(
    Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
    std::vector<std::pair<Graph, Graph>> mats_keys2elems, TagBase const* tagbase,
    LOs same_ents2old_ents, LOs same_ents2new_ents, bool conserves_density) {
  auto dim = old_mesh->dim();
  auto ncomps = tagbase->ncomps();
  auto velocity_name = tagbase->name();
  auto momentum_name = opts.velocity_momentum_map.find(velocity_name).second
  auto density_name = opts.velocity_density_map.find(velocity_name).second
  auto error_name = momentum_name + "_error";
  auto old_elem_densities = old_mesh->get_array<Real>(dim, density_name);
  auto new_elem_densities = new_mesh->get_array<Real>(dim, density_name);
  auto old_vert_velocities = old_mesh->get_array<Real>(VERT, velocity_name);
  auto new_vert_velocities = new_mesh->get_array<Real>(VERT, velocity_name);
  auto old_elem_velocities = average_field(old_mesh, dim, ncomps, old_vert_velocities);
  auto new_elem_velocities = average_field(new_mesh, dim, ncomps, new_vert_velocities);
  auto old_elem_momenta = multiply_each(old_elem_velocities, old_elem_densities);
  auto new_elem_momenta = multiply_each(new_elem_velocities, new_elem_densities);
  track_cavs_integral_error(old_mesh, new_mesh, mats_keys2elems,
      old_elem_momenta, new_elem_momenta, error_name, conserves_density);
}

/* we separate sub-cavities even further than by material: we also separate
   the elements touching the boundary from those that don't.
   this ensures that size error, which can only be created and repaired at
   the boundary, remains at the boundary.
   it should also be beneficial to do this for density errors, because
   a boundary collapse creates a correspondingly large integral error
   for each size error, and so those integral errors will be corrected
   in the same elements in which the size errors are corrected.
 */
static LOs negate_boundary_elem_class_ids(Mesh* mesh) {
  auto vert_class_dims = mesh->get_array<I8>(VERT, "class_dim");
  auto verts_are_int = each_eq_to(vert_class_dims, I8(mesh->dim()));
  auto verts_are_bdry = invert_marks(verts_are_int);
  auto elems_are_bdry = mark_up(mesh, VERT, mesh->dim(), verts_are_bdry);
  LOs elem_class_ids;
  if (mesh->has_tag<LO>(mesh->dim(), "class_id")) {
    elem_class_ids = mesh->get_tag<LO>(mesh->dim(), "class_id");
  } else {
    elem_class_ids = LOs(mesh->nelems(), 1);
  }
  auto out = deep_copy(elem_class_ids);
  auto f = LAMBDA(LO e) {
    out[e] = (elems_are_bdry[e]) ? (-elem_class_ids[e]) : (elem_class_ids[e]);
  };
  parallel_for(mesh->nelems(), f);
  return out;
}

void track_cavs_all_errors(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
    Int key_dim, LOs keys2kds, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents,
    bool conserves_size, bool conserves_mass, bool conserves_momentum) {
  auto dim = new_mesh->dim();
  auto kds2old_elems = old_mesh->ask_up(key_dim, dim);
  auto keys2old_elems = unmap_graph(keys2kds, kds2old_elems);
  auto keys2new_elems = Graph(keys2prods, prods2new_ents);
  std::vector<std::pair<Graph, Graph>> mats_keys2elems;
  auto old_class_ids = negate_boundary_elem_class_ids(old_mesh);
  auto new_class_ids = negate_boundary_elem_class_ids(new_mesh);
  mats_keys2elems = separate_cavities(
      keys2old_elems, old_class_ids, keys2new_elems, new_class_ids);
  if (opts.should_conserve_size) {
    track_cavs_size_error(old_mesh, new_mesh, mats_keys2elems,
        same_ents2old_ents, same_ents2new_ents, can_create_errors);
  }
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_conserve(old_mesh, opts, dim, tagbase)) {
      track_cavs_density_error(old_mesh, opts, new_mesh,
          mats_keys2elems, tagbase,
          same_ents2old_ents, same_ents2new_ents, can_create_errors);
    }
  }
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tagbase = old_mesh->get_tag(VERT, i);
    if (is_momentum_velocity(old_mesh, opts, VERT, tagbase)) {
      track_cavs_momentum_error(old_mesh, opts, new_mesh,
          mats_keys2elems, tagbase,
          same_ents2old_ents, same_ents2new_ents, can_create_errors);
    }
  }
}

void transfer_conserve_refine(Mesh* old_mesh, XferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  if (!should_conserve_any(old_mesh, opts)) return;
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_conserve(old_mesh, opts, dim, tagbase)) {
      /* just inherit the density field */
      transfer_inherit_refine(old_mesh, new_mesh, keys2edges, old_mesh->dim(),
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase);
    }
  }
  bool conserves_size = true;
  bool conserves_mass = true;
  bool conserves_momentum = false;
  track_cavs_all_errors(old_mesh, opts, new_mesh, EDGE, keys2edges, keys2prods,
      prods2new_ents, same_ents2old_ents, same_ents2new_ents,
      conserves_size, conserves_mass, conserves_momentum);
}

/* intersection-based transfer of density fields.
   note that this exactly conserves integrals on the interior,
   but we purposefully omit the modifications that would
   ensure conservation for boundary collapses that change the
   domain. instead, those errors will be tracked by the code above. */
template <Int dim>
static void transfer_subcavs_density_dim(Mesh* old_mesh, Mesh* new_mesh,
    TagBase const* tagbase, Graph keys2old_elems, Graph keys2new_elems,
    Write<Real> new_data_w) {
  auto ncomps = tagbase->ncomps();
  auto old_tag = to<Real>(tagbase);
  auto old_data = old_tag->array();
  auto old_ev2v = old_mesh->ask_elem_verts();
  auto old_coords = old_mesh->coords();
  auto new_ev2v = new_mesh->ask_elem_verts();
  auto new_coords = new_mesh->coords();
  auto new_sizes = new_mesh->ask_sizes();
  auto nkeys = keys2old_elems.nnodes();
  auto f = LAMBDA(LO key) {
    for (auto kne = keys2new_elems.a2ab[key];
         kne < keys2new_elems.a2ab[key + 1]; ++kne) {
      auto new_elem = keys2new_elems.ab2b[kne];
      for (Int comp = 0; comp < ncomps; ++comp) {
        new_data_w[new_elem * ncomps + comp] = 0;
      }
      auto new_verts = gather_verts<dim + 1>(new_ev2v, new_elem);
      auto new_points =
          gather_vectors<dim + 1, dim>(new_coords, new_verts);
      for (auto koe = keys2old_elems.a2ab[key];
           koe < keys2old_elems.a2ab[key + 1]; ++koe) {
        auto old_elem = keys2old_elems.ab2b[koe];
        auto old_verts = gather_verts<dim + 1>(old_ev2v, old_elem);
        auto old_points = gather_vectors<dim + 1, dim>(old_coords, old_verts);
        r3d::Polytope<dim> intersection;
        r3d::intersect_simplices(
            intersection, to_r3d(target_points), to_r3d(donor_points));
        auto intersection_size = r3d::measure(intersection);
        for (Int comp = 0; comp < ncomps; ++comp) {
          new_data_w[new_elem * ncomps + comp] +=
            intersection_size * old_data[donor_elem * ncomps + comp];
        }
      }
      for (Int comp = 0; comp < ncomps; ++comp) {
        new_data_w[new_elems * ncomps + comp] /= new_sizes[new_elem]
      }
    }
  };
  parallel_for(nkeys, f);
}

static void transfer_cavs_density(Mesh* old_mesh, Mesh* new_mesh,
    std::vector<std::pair<Graph, Graph>> mats_keys2elems,
    LOs same_ents2old_ents, LOs same_ents2new_ents, TagBase const* tagbase) {
  auto dim = old_mesh->dim();
  auto nnew_elems = new_mesh->nelems();
  auto ncomps = tagbase->ncomps();
  auto new_data_w = Write<Real>(nnew_elems * ncomps);
  for (auto pair : mats_keys2elems) {
    auto keys2old_elems = pair.first;
    auto keys2new_elems = pair.second;
    if (dim == 3) {
      transfer_conserve_dim<3>(old_mesh, new_mesh, tagbase, keys2old_elems,
          keys2new_elems, new_data_w);
    } else if (dim == 2) {
      transfer_conserve_dim<2>(old_mesh, new_mesh, tagbase, keys2old_elems,
          keys2new_elems, new_data_w);
    }
  }
  transfer_common2(old_mesh, new_mesh, dim, same_ents2old_ents,
      same_ents2new_ents, tagbase, new_data_w);
}

void transfer_cavs_all_densities(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
    Int key_dim, LOs keys2kds, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  auto dim = new_mesh->dim();
  /* TODO: consolidate this graph creation code with the one in
     track_cavs_all_errors() */
  auto kds2old_elems = old_mesh->ask_up(key_dim, dim);
  auto keys2old_elems = unmap_graph(keys2kds, kds2old_elems);
  auto keys2new_elems = Graph(keys2prods, prods2new_ents);
  std::vector<std::pair<Graph, Graph>> mats_keys2elems;
  if (old_mesh->has_tag(VERT, "class_id")) {
    auto old_class_ids = old_mesh->get_array<LO>(VERT, "class_id");
    auto new_class_ids = new_mesh->get_array<LO>(VERT, "class_id");
    mats_keys2elems = separate_cavities(
        keys2old_elems, old_class_ids, keys2new_elems, new_class_ids);
  } else {
    mats_keys2elems = {{keys2old_elems, keys2new_elems}};
  }
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_conserve(old_mesh, opts, dim, tagbase)) {
      transfer_cavs_density(old_mesh, new_mesh,
          mats_keys2elems,
          same_ents2old_ents, same_ents2new_ents, tagbase);
    }
  }
}

void transfer_conserve_swap(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
    LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  transfer_cavs_all_densities(old_mesh, opts, new_mesh, key_dim,
      keys2kds, keys2prods, prods2new_ents,
      same_ents2old_ents, same_ents2new_ents);
  bool conserves_size = true;
  bool conserves_mass = true;
  bool conserves_momentum = false;
  track_cavs_all_errors(old_mesh, opts, new_mesh, EDGE, keys2edges, keys2prods,
      prods2new_ents, same_ents2old_ents, same_ents2new_ents,
      conserves_size, conserves_mass, conserves_momentum);
}

/* conservative diffusion:
   1) diffusion is not allowed across classification boundaries.
      this preserves conservation on a per-object basis
   2) we are careful with the mathematics of the diffusion algorithm,
      such that it should preserve the total error, and hence
      conservation
 */
static Graph get_elem_diffusion_graph(Mesh* mesh) {
  auto dim = mesh->dim();
  OMEGA_H_CHECK(mesh->owners_have_all_upward(dim - 1));
  auto side_class_dims = mesh->get_array<I8>(dim - 1, "class_dim");
  auto sides_are_int = each_eq_to(side_class_dims, I8(dim));
  auto sides_are_bdry = invert_marks(sides_are_int);
  auto elems2sides = mesh->ask_down(dim, dim - 1);
  auto sides2elems = mesh->ask_up(dim - 1, dim);
  return elements_across_sides(dim, elems2sides, sides2elems, sides_are_bdry);
}

static Reals diffuse_elem_error_once(Mesh* mesh, Graph g, Reals x, Int ncomps) {
  auto out = deep_copy(x);
  auto max_deg = mesh->dim() + 1;
  auto f = OMEGA_H_LAMBDA(LO e) {
    for (auto ee = g.a2ab[e]; ee < g.a2ab[e + 1]; ++ee) {
      auto oe = g.ab2b[ee];
      for (Int c = 0; c < ncomps; ++c) {
        out[e * ncomps + c] += (x[oe * ncomps + c] - x[e * ncomps + c]) / max_deg;
      }
    }
  };
  parallel_for(mesh->nelems(), f);
  return mesh->sync_array(mesh->dim(), out, ncomps);
}

static void diffuse_elem_error_tag_once(Mesh* mesh, Graph g, std::string const& error_name) {
  auto array = mesh->get_array<Real>(dim, error_name);
  auto ncomps = tagbase->ncomps();
  array = diffuse_elem_error_once(mesh, g, array, ncomps);
  mesh->set_tag<Real>(dim, error_name, array);
}

static void diffuse_all_errors(Mesh* mesh, XferOpts const& opts) {
  auto g = get_elem_diffusion_graph(mesh);
  auto dim = mesh->dim();
  auto niters = 8;
  for (Int iter = 0; iter < niters; ++iter) {
    for (Int tagi = 0; tagi < old_mesh->ntags(dim); ++tagi) {
      auto tagbase = mesh->get_tag(dim, tagi);
      if (should_conserve(mesh, opts, dim, tagbase)) {
        auto integral_name = opts.integral_map.find(tagbase->name()).second;
        auto error_name = integral_name + "_error";
        diffuse_elem_error_tag_once(mesh, g, error_name);
      }
    }
    for (Int tagi = 0; tagi < old_mesh->ntags(VERT); ++tagi) {
      auto tagbase = mesh->get_tag(VERT, tagi);
      if (should_conserve(mesh, opts, VERT, tagbase)) {
        auto integral_name = opts.integral_map.find(tagbase->name()).second;
        auto error_name = integral_name + "_error";
        diffuse_elem_error_tag_once(mesh, g, error_name);
      }
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

static Read<I8> get_comps_are_fixed(Mesh* mesh) {
  if (mesh->has_tag(VERT, "momentum_velocity_fixed")) {
    return mesh->get_array<I8>(VERT, "momentum_velocity_fixed");
  } else {
    return Read<I8>(mesh->nverts(), I8(0));
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

}  // end namespace Omega_h
