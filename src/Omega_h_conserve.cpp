#include "Omega_h_conserve.hpp"

#include "Omega_h_map.hpp"

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
  auto conserves_mass = true;
  auto conserves_mass = true;
  track_cavs_all_errors(old_mesh, opts, new_mesh,
      EDGE, keys2edges, keys2prods, prods2new_ents,
      same_ents2old_ents, same_ents2new_ents, can_create_errors);
}

void transfer_conserve_swap(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
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
  if (opts.should_conserve_size) {
    track_cav_size_error(old_mesh, new_mesh, mats_keys2elems,
        same_ents2old_ents, same_ents2new_ents);
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
