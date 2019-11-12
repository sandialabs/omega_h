#include "Omega_h_transfer.hpp"

#include "Omega_h_affine.hpp"
#include "Omega_h_conserve.hpp"
#include "Omega_h_fit.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_shape.hpp"

namespace Omega_h {

bool is_transfer_required(
    TransferOpts const& opts, std::string const& name, Omega_h_Transfer type) {
  if (!opts.type_map.count(name)) return false;
  return opts.type_map.find(name)->second == type;
}

bool should_inherit(
    Mesh* mesh, TransferOpts const& opts, Int, TagBase const* tag) {
  auto& name = tag->name();
  if (!(is_transfer_required(opts, name, OMEGA_H_INHERIT) ||
          name == "class_id" || name == "class_dim" ||
          name == "momentum_velocity_fixed")) {
    return false;
  }
  for (Int i = 0; i <= mesh->dim(); ++i) {
    if (!mesh->has_tag(i, name)) return false;
    if (mesh->get_tagbase(i, name)->type() != tag->type()) return false;
    if (mesh->get_tagbase(i, name)->ncomps() != tag->ncomps()) return false;
  }
  return true;
}

bool should_interpolate(
    Mesh*, TransferOpts const& opts, Int dim, TagBase const* tag) {
  auto& name = tag->name();
  if (!(is_transfer_required(opts, name, OMEGA_H_LINEAR_INTERP) ||
          is_transfer_required(opts, name, OMEGA_H_MOMENTUM_VELOCITY) ||
          name == "coordinates" || name == "warp")) {
    return false;
  }
  return dim == VERT && tag->type() == OMEGA_H_REAL;
}

bool should_fit(
    Mesh* mesh, TransferOpts const& opts, Int dim, TagBase const* tag) {
  auto& name = tag->name();
  if (!is_transfer_required(opts, name, OMEGA_H_POINTWISE)) {
    return false;
  }
  return dim == mesh->dim() && tag->type() == OMEGA_H_REAL;
}

bool is_density(
    Mesh* mesh, TransferOpts const& opts, Int dim, TagBase const* tag) {
  auto& name = tag->name();
  if (!is_transfer_required(opts, name, OMEGA_H_DENSITY)) {
    return false;
  }
  return dim == mesh->dim() && tag->type() == OMEGA_H_REAL;
}

bool should_conserve(
    Mesh* mesh, TransferOpts const& opts, Int dim, TagBase const* tag) {
  auto& name = tag->name();
  if (!is_transfer_required(opts, name, OMEGA_H_CONSERVE)) {
    return false;
  }
  return dim == mesh->dim() && tag->type() == OMEGA_H_REAL;
}

bool has_densities_or_conserved(Mesh* mesh, TransferOpts const& opts) {
  auto dim = mesh->dim();
  for (Int i = 0; i < mesh->ntags(dim); ++i) {
    auto tag = mesh->get_tag(dim, i);
    if (should_conserve(mesh, opts, dim, tag)) {
      return true;
    }
    if (is_density(mesh, opts, dim, tag)) {
      return true;
    }
  }
  return false;
}

bool should_conserve_any(Mesh* mesh, TransferOpts const& opts) {
  auto dim = mesh->dim();
  for (Int i = 0; i < mesh->ntags(dim); ++i) {
    if (should_conserve(mesh, opts, dim, mesh->get_tag(dim, i))) {
      return true;
    }
  }
  return false;
}

bool is_metric(
    Mesh* mesh, TransferOpts const& opts, Int dim, TagBase const* tag) {
  auto& name = tag->name();
  if (!(is_transfer_required(opts, name, OMEGA_H_METRIC) || name == "metric" ||
          name == "target_metric")) {
    return false;
  }
  return dim == VERT && tag->type() == OMEGA_H_REAL &&
         (tag->ncomps() == symm_ncomps(mesh->dim()) || tag->ncomps() == 1);
}

bool is_momentum_velocity(
    Mesh* mesh, TransferOpts const& opts, Int dim, TagBase const* tag) {
  auto& name = tag->name();
  if (!is_transfer_required(opts, name, OMEGA_H_MOMENTUM_VELOCITY)) {
    return false;
  }
  return dim == VERT && tag->type() == OMEGA_H_REAL &&
         tag->ncomps() == mesh->dim();
}

bool has_momentum_velocity(Mesh* mesh, TransferOpts const& opts) {
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    if (is_momentum_velocity(mesh, opts, VERT, mesh->get_tag(VERT, i))) {
      return true;
    }
  }
  return false;
}

static bool should_transfer_no_products(
    Mesh* mesh, TransferOpts const& opts, Int dim, TagBase const* tag) {
  return should_inherit(mesh, opts, dim, tag) ||
         should_interpolate(mesh, opts, dim, tag) ||
         is_metric(mesh, opts, dim, tag) ||
         should_conserve(mesh, opts, dim, tag) ||
         is_density(mesh, opts, dim, tag) ||
         is_momentum_velocity(mesh, opts, dim, tag);
}

static bool should_transfer_copy(
    Mesh* mesh, TransferOpts const& opts, Int dim, TagBase const* tag) {
  return should_transfer_no_products(mesh, opts, dim, tag) ||
         tag->name() == "global";
}

static bool should_transfer_density(
    Mesh* mesh, TransferOpts const& opts, Int dim, TagBase const* tag) {
  return should_conserve(mesh, opts, dim, tag) ||
         is_density(mesh, opts, dim, tag);
}

template <typename T>
void transfer_common3(
    Mesh* new_mesh, Int ent_dim, TagBase const* tagbase, Write<T> new_data) {
  auto const& name = tagbase->name();
  auto ncomps = tagbase->ncomps();
  new_mesh->add_tag(ent_dim, name, ncomps, Read<T>(new_data), true);
}

template <typename T>
void transfer_common2(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs same_ents2old_ents, LOs same_ents2new_ents, TagBase const* tagbase,
    Write<T> new_data) {
  auto const& name = tagbase->name();
  auto ncomps = tagbase->ncomps();
  auto old_data = old_mesh->get_array<T>(ent_dim, name);
  auto same_data = read(unmap(same_ents2old_ents, old_data, ncomps));
  map_into(same_data, same_ents2new_ents, new_data, ncomps);
  transfer_common3(new_mesh, ent_dim, tagbase, new_data);
}

template <typename T>
void transfer_common(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    LOs same_ents2old_ents, LOs same_ents2new_ents, LOs prods2new_ents,
    TagBase const* tagbase, Read<T> prod_data) {
  auto nnew_ents = new_mesh->nents(ent_dim);
  auto ncomps = tagbase->ncomps();
  auto new_data = Write<T>(nnew_ents * ncomps);
  map_into(prod_data, prods2new_ents, new_data, ncomps);
  transfer_common2(old_mesh, new_mesh, ent_dim, same_ents2old_ents,
      same_ents2new_ents, tagbase, new_data);
}

static void transfer_linear_interp(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2midverts, LOs same_verts2old_verts,
    LOs same_verts2new_verts) {
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tagbase = old_mesh->get_tag(VERT, i);
    if (should_interpolate(old_mesh, opts, VERT, tagbase)) {
      auto ncomps = tagbase->ncomps();
      auto old_data = old_mesh->get_array<Real>(VERT, tagbase->name());
      auto prod_data =
          average_field(old_mesh, EDGE, keys2edges, ncomps, old_data);
      transfer_common(old_mesh, new_mesh, VERT, same_verts2old_verts,
          same_verts2new_verts, keys2midverts, tagbase, prod_data);
    }
  }
}

static void transfer_metric(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2midverts, LOs same_verts2old_verts,
    LOs same_verts2new_verts) {
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tagbase = old_mesh->get_tag(VERT, i);
    if (is_metric(old_mesh, opts, VERT, tagbase)) {
      auto old_data = old_mesh->get_array<Real>(VERT, tagbase->name());
      auto prod_data = get_mident_metrics(old_mesh, EDGE, keys2edges, old_data);
      transfer_common(old_mesh, new_mesh, VERT, same_verts2old_verts,
          same_verts2new_verts, keys2midverts, tagbase, prod_data);
    }
  }
}

template <typename T>
void transfer_inherit_refine(Mesh* old_mesh, Mesh* new_mesh, LOs keys2edges,
    Int prod_dim, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents, std::string const& name) {
  auto old_tag = old_mesh->get_tag<T>(prod_dim, name);
  auto ncomps = old_tag->ncomps();
  auto nprods = keys2prods.last();
  auto prod_data = Write<T>(nprods * ncomps);
  auto nkeys = keys2edges.size();
  /* transfer pairs */
  if (prod_dim > VERT) {
    auto dom_dim = prod_dim;
    auto dom_data = old_mesh->get_array<T>(dom_dim, name);
    auto edges2doms = old_mesh->ask_graph(EDGE, dom_dim);
    auto edges2edge_doms = edges2doms.a2ab;
    auto edge_doms2doms = edges2doms.ab2b;
    auto f = OMEGA_H_LAMBDA(LO key) {
      auto edge = keys2edges[key];
      auto prod = keys2prods[key];
      for (auto edge_dom = edges2edge_doms[edge];
           edge_dom < edges2edge_doms[edge + 1]; ++edge_dom) {
        auto dom = edge_doms2doms[edge_dom];
        for (Int pair = 0; pair < 2; ++pair) {
          for (Int comp = 0; comp < ncomps; ++comp) {
            prod_data[prod * ncomps + comp] = dom_data[dom * ncomps + comp];
          }
          ++prod;
        }
      }
    };
    parallel_for(nkeys, f, "transfer_inherit_refine(pairs)");
  }
  if (prod_dim < old_mesh->dim()) {
    auto dom_dim = prod_dim + 1;
    auto dom_data = old_mesh->get_array<T>(dom_dim, name);
    auto edges2doms = old_mesh->ask_graph(EDGE, dom_dim);
    auto edges2edge_doms = edges2doms.a2ab;
    auto edge_doms2doms = edges2doms.ab2b;
    auto f = OMEGA_H_LAMBDA(LO key) {
      auto edge = keys2edges[key];
      auto ndoms = edges2edge_doms[edge + 1] - edges2edge_doms[edge];
      auto prod = keys2prods[key + 1] - ndoms;
      for (auto edge_dom = edges2edge_doms[edge];
           edge_dom < edges2edge_doms[edge + 1]; ++edge_dom) {
        auto dom = edge_doms2doms[edge_dom];
        for (Int comp = 0; comp < ncomps; ++comp) {
          prod_data[prod * ncomps + comp] = dom_data[dom * ncomps + comp];
        }
        ++prod;
      }
    };
    parallel_for(nkeys, f, "transfer_inherit_refine(cuts)");
  }
  transfer_common(old_mesh, new_mesh, prod_dim, same_ents2old_ents,
      same_ents2new_ents, prods2new_ents, old_tag, Read<T>(prod_data));
}

void transfer_inherit_refine(Mesh* old_mesh, Mesh* new_mesh, LOs keys2edges,
    Int prod_dim, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents, TagBase const* tagbase) {
  switch (tagbase->type()) {
    case OMEGA_H_I8:
      transfer_inherit_refine<I8>(old_mesh, new_mesh, keys2edges, prod_dim,
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase->name());
      break;
    case OMEGA_H_I32:
      transfer_inherit_refine<I32>(old_mesh, new_mesh, keys2edges, prod_dim,
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase->name());
      break;
    case OMEGA_H_I64:
      transfer_inherit_refine<I64>(old_mesh, new_mesh, keys2edges, prod_dim,
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase->name());
      break;
    case OMEGA_H_F64:
      transfer_inherit_refine<Real>(old_mesh, new_mesh, keys2edges, prod_dim,
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase->name());
      break;
  }
}

static void transfer_inherit_refine(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, Int prod_dim, LOs keys2prods,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents) {
  for (Int i = 0; i < old_mesh->ntags(prod_dim); ++i) {
    auto tagbase = old_mesh->get_tag(prod_dim, i);
    if (should_inherit(old_mesh, opts, prod_dim, tagbase)) {
      transfer_inherit_refine(old_mesh, new_mesh, keys2edges, prod_dim,
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase);
    }
  }
}

static void transfer_density_refine(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_transfer_density(old_mesh, opts, dim, tagbase)) {
      /* just inherit the density field */
      transfer_inherit_refine(old_mesh, new_mesh, keys2edges, old_mesh->dim(),
          keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
          tagbase);
    }
  }
}

static void transfer_pointwise_refine(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_fit(old_mesh, opts, dim, tagbase)) {
      transfer_inherit_refine(old_mesh, new_mesh, keys2edges, dim, keys2prods,
          prods2new_ents, same_ents2old_ents, same_ents2new_ents, tagbase);
    }
  }
}

void transfer_length(Mesh* old_mesh, Mesh* new_mesh, LOs same_ents2old_ents,
    LOs same_ents2new_ents, LOs prods2new_ents) {
  for (Int i = 0; i < old_mesh->ntags(EDGE); ++i) {
    auto tagbase = old_mesh->get_tag(EDGE, i);
    if (tagbase->name() == "length" && tagbase->type() == OMEGA_H_REAL &&
        tagbase->ncomps() == 1) {
      auto prod_data = measure_edges_metric(new_mesh, prods2new_ents);
      transfer_common(old_mesh, new_mesh, EDGE, same_ents2old_ents,
          same_ents2new_ents, prods2new_ents, tagbase, prod_data);
    }
  }
}

void transfer_quality(Mesh* old_mesh, Mesh* new_mesh, LOs same_ents2old_ents,
    LOs same_ents2new_ents, LOs prods2new_ents) {
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (tagbase->name() == "quality" && tagbase->type() == OMEGA_H_REAL &&
        tagbase->ncomps() == 1) {
      auto prod_data = measure_qualities(new_mesh, prods2new_ents);
      transfer_common(old_mesh, new_mesh, dim, same_ents2old_ents,
          same_ents2new_ents, prods2new_ents, tagbase, prod_data);
    }
  }
}

void transfer_size(Mesh* old_mesh, Mesh* new_mesh, LOs same_ents2old_ents,
    LOs same_ents2new_ents, LOs prods2new_ents) {
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (tagbase->name() == "size" && tagbase->type() == OMEGA_H_REAL &&
        tagbase->ncomps() == 1) {
      auto prod_data = measure_elements_real(new_mesh, prods2new_ents);
      transfer_common(old_mesh, new_mesh, dim, same_ents2old_ents,
          same_ents2new_ents, prods2new_ents, tagbase, prod_data);
    }
  }
}

static void transfer_face_flux(Mesh* old_mesh, Mesh* new_mesh,
    LOs same_ents2old_ents, LOs same_ents2new_ents, LOs prods2new_ents) {
  for (Int i = 0; i < old_mesh->ntags(FACE); ++i) {
    TagBase const* tagbase = old_mesh->get_tag(FACE, i);
    if (tagbase->name() == "magnetic face flux") {
      Read<Real> old_data = old_mesh->get_array<Real>(FACE, tagbase->name());
      Read<Real> prod_data(prods2new_ents.size(), 0);
      transfer_common(old_mesh, new_mesh, FACE, same_ents2old_ents,
          same_ents2new_ents, prods2new_ents, tagbase, prod_data);
    }
  }
}

void transfer_refine(Mesh* old_mesh, TransferOpts const& opts, Mesh* new_mesh,
    LOs keys2edges, LOs keys2midverts, Int prod_dim, LOs keys2prods,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents) {
  begin_code("transfer_refine");
  transfer_inherit_refine(old_mesh, opts, new_mesh, keys2edges, prod_dim,
      keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents);
  if (prod_dim == VERT) {
    transfer_linear_interp(old_mesh, opts, new_mesh, keys2edges, keys2midverts,
        same_ents2old_ents, same_ents2new_ents);
    transfer_metric(old_mesh, opts, new_mesh, keys2edges, keys2midverts,
        same_ents2old_ents, same_ents2new_ents);
  }
  if (prod_dim == EDGE) {
    transfer_length(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
  } else if (prod_dim == FACE) {
    transfer_face_flux(old_mesh, new_mesh, same_ents2old_ents,
        same_ents2new_ents, prods2new_ents);
  }
  if (prod_dim == old_mesh->dim()) {
    transfer_size(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
    transfer_quality(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
    transfer_density_refine(old_mesh, opts, new_mesh, keys2edges, keys2prods,
        prods2new_ents, same_ents2old_ents, same_ents2new_ents);
    transfer_conserve_refine(old_mesh, opts, new_mesh, keys2edges, keys2prods,
        prods2new_ents, same_ents2old_ents, same_ents2new_ents);
    transfer_pointwise_refine(old_mesh, opts, new_mesh, keys2edges, keys2prods,
        prods2new_ents, same_ents2old_ents, same_ents2new_ents);
  }
  if (opts.user_xfer) {
    opts.user_xfer->refine(*old_mesh, *new_mesh, keys2edges, keys2midverts,
        prod_dim, keys2prods, prods2new_ents, same_ents2old_ents,
        same_ents2new_ents);
  }
  end_code();
}

template <typename T>
static void transfer_inherit_coarsen_tmpl(Mesh* old_mesh, Mesh* new_mesh,
    Adj keys2doms, Int prod_dim, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents, TagBase const* tagbase) {
  auto name = tagbase->name();
  auto old_tag = as<T>(tagbase);
  auto ncomps = old_tag->ncomps();
  auto dom_data = old_tag->array();
  auto key_doms2doms = keys2doms.ab2b;
  auto prod_data = read(unmap(key_doms2doms, dom_data, ncomps));
  transfer_common(old_mesh, new_mesh, prod_dim, same_ents2old_ents,
      same_ents2new_ents, prods2new_ents, old_tag, prod_data);
}

static void transfer_inherit_coarsen(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, Adj keys2doms, Int prod_dim, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  for (Int i = 0; i < old_mesh->ntags(prod_dim); ++i) {
    auto tagbase = old_mesh->get_tag(prod_dim, i);
    if (should_inherit(old_mesh, opts, prod_dim, tagbase)) {
      switch (tagbase->type()) {
        case OMEGA_H_I8:
          transfer_inherit_coarsen_tmpl<I8>(old_mesh, new_mesh, keys2doms,
              prod_dim, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
              tagbase);
          break;
        case OMEGA_H_I32:
          transfer_inherit_coarsen_tmpl<I32>(old_mesh, new_mesh, keys2doms,
              prod_dim, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
              tagbase);
          break;
        case OMEGA_H_I64:
          transfer_inherit_coarsen_tmpl<I64>(old_mesh, new_mesh, keys2doms,
              prod_dim, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
              tagbase);
          break;
        case OMEGA_H_F64:
          transfer_inherit_coarsen_tmpl<Real>(old_mesh, new_mesh, keys2doms,
              prod_dim, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
              tagbase);
          break;
      }
    }
  }
}

template <typename T>
static void transfer_no_products_tmpl(Mesh* old_mesh, Mesh* new_mesh,
    Int prod_dim, LOs same_ents2old_ents, LOs same_ents2new_ents,
    TagBase const* tagbase) {
  auto old_tag = as<T>(tagbase);
  auto prods2new_ents = LOs({});
  auto prod_data = Read<T>({});
  transfer_common(old_mesh, new_mesh, prod_dim, same_ents2old_ents,
      same_ents2new_ents, prods2new_ents, old_tag, prod_data);
}

static void transfer_no_products(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, Int prod_dim, LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  for (Int i = 0; i < old_mesh->ntags(prod_dim); ++i) {
    auto tagbase = old_mesh->get_tag(prod_dim, i);
    if (should_transfer_no_products(old_mesh, opts, prod_dim, tagbase)) {
      switch (tagbase->type()) {
        case OMEGA_H_I8:
          transfer_no_products_tmpl<I8>(old_mesh, new_mesh, prod_dim,
              same_ents2old_ents, same_ents2new_ents, tagbase);
          break;
        case OMEGA_H_I32:
          transfer_no_products_tmpl<I32>(old_mesh, new_mesh, prod_dim,
              same_ents2old_ents, same_ents2new_ents, tagbase);
          break;
        case OMEGA_H_I64:
          transfer_no_products_tmpl<I64>(old_mesh, new_mesh, prod_dim,
              same_ents2old_ents, same_ents2new_ents, tagbase);
          break;
        case OMEGA_H_F64:
          transfer_no_products_tmpl<Real>(old_mesh, new_mesh, prod_dim,
              same_ents2old_ents, same_ents2new_ents, tagbase);
          break;
      }
    }
  }
}

template <Int dim>
void transfer_pointwise_tmpl(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_elems, LOs same_elems2old_elems,
    LOs same_elems2new_elems, TagBase const* tagbase) {
  auto name = tagbase->name();
  auto old_tag = as<Real>(tagbase);
  auto ncomps = old_tag->ncomps();
  auto old_data = old_tag->array();
  auto kds2elems = old_mesh->ask_up(key_dim, dim);
  auto kds2kd_elems = kds2elems.a2ab;
  auto kd_elems2elems = kds2elems.ab2b;
  auto old_elem_verts2verts = old_mesh->ask_verts_of(dim);
  auto old_coords = old_mesh->coords();
  auto new_elem_verts2verts = new_mesh->ask_verts_of(dim);
  auto new_coords = new_mesh->coords();
  auto nkeys = keys2kds.size();
  auto nprods = keys2prods.last();
  auto prod_data_w = Write<Real>(nprods * ncomps);
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto kd = keys2kds[key];
    for (auto prod = keys2prods[key]; prod < keys2prods[key + 1]; ++prod) {
      auto new_elem = prods2new_elems[prod];
      auto new_v = gather_verts<dim + 1>(new_elem_verts2verts, new_elem);
      auto new_vp = gather_vectors<dim + 1, dim>(new_coords, new_v);
      auto new_ip = average(new_vp);
      LO best_old_elem = -1;
      Real best_distance = 0.0;
      for (auto kde = kds2kd_elems[kd]; kde < kds2kd_elems[kd + 1]; ++kde) {
        auto old_elem = kd_elems2elems[kde];
        auto old_v = gather_verts<dim + 1>(old_elem_verts2verts, old_elem);
        auto old_vp = gather_vectors<dim + 1, dim>(old_coords, old_v);
        auto a = simplex_affine(old_vp);
        auto ia = invert(a);
        auto xi = ia * new_ip;
        auto bc = form_barycentric(xi);
        auto distance = -reduce(bc, minimum<Real>());
        if (best_old_elem == -1 || distance < best_distance) {
          best_old_elem = old_elem;
          best_distance = distance;
        }
      }
      for (Int comp = 0; comp < ncomps; ++comp) {
        prod_data_w[prod * ncomps + comp] =
            old_data[best_old_elem * ncomps + comp];
      }
    }
  };
  parallel_for(nkeys, f, "transfer_pointwise");
  auto prod_data = Reals(prod_data_w);
  transfer_common(old_mesh, new_mesh, dim, same_elems2old_elems,
      same_elems2new_elems, prods2new_elems, old_tag, prod_data);
}

void transfer_pointwise(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, Int key_dim, LOs keys2kds, LOs keys2prods,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents) {
  auto dim = new_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_fit(old_mesh, opts, dim, tagbase)) {
      if (dim == 3) {
        transfer_pointwise_tmpl<3>(old_mesh, new_mesh, key_dim, keys2kds,
            keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
            tagbase);
      } else if (dim == 2) {
        transfer_pointwise_tmpl<2>(old_mesh, new_mesh, key_dim, keys2kds,
            keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents,
            tagbase);
      }
    }
  }
}

void transfer_coarsen(Mesh* old_mesh, TransferOpts const& opts, Mesh* new_mesh,
    LOs keys2verts, Adj keys2doms, Int prod_dim, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  begin_code("transfer_coarsen");
  if (prod_dim == VERT) {
    transfer_no_products(old_mesh, opts, new_mesh, prod_dim, same_ents2old_ents,
        same_ents2new_ents);
  } else {
    transfer_inherit_coarsen(old_mesh, opts, new_mesh, keys2doms, prod_dim,
        prods2new_ents, same_ents2old_ents, same_ents2new_ents);
  }
  if (prod_dim == EDGE) {
    transfer_length(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
  } else if (prod_dim == FACE) {
    transfer_face_flux(old_mesh, new_mesh, same_ents2old_ents,
        same_ents2new_ents, prods2new_ents);
  }
  if (prod_dim == old_mesh->dim()) {
    transfer_size(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
    transfer_quality(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
    transfer_pointwise(old_mesh, opts, new_mesh, VERT, keys2verts,
        keys2doms.a2ab, prods2new_ents, same_ents2old_ents, same_ents2new_ents);
    transfer_densities_and_conserve_coarsen(old_mesh, opts, new_mesh,
        keys2verts, keys2doms, prods2new_ents, same_ents2old_ents,
        same_ents2new_ents);
  }
  if (opts.user_xfer) {
    opts.user_xfer->coarsen(*old_mesh, *new_mesh, keys2verts, keys2doms,
        prod_dim, prods2new_ents, same_ents2old_ents, same_ents2new_ents);
  }
  end_code();
}

template <typename T>
static void transfer_copy_tmpl(
    Mesh* new_mesh, Int prod_dim, TagBase const* tagbase) {
  auto old_tag = as<T>(tagbase);
  auto const& name = old_tag->name();
  auto ncomps = old_tag->ncomps();
  auto old_data = old_tag->array();
  new_mesh->add_tag(prod_dim, name, ncomps, old_data, true);
}

void transfer_copy(
    Mesh* old_mesh, TransferOpts const& opts, Mesh* new_mesh, Int prod_dim) {
  for (Int i = 0; i < old_mesh->ntags(prod_dim); ++i) {
    auto tagbase = old_mesh->get_tag(prod_dim, i);
    if (should_transfer_copy(old_mesh, opts, prod_dim, tagbase)) {
      switch (tagbase->type()) {
        case OMEGA_H_I8:
          transfer_copy_tmpl<I8>(new_mesh, prod_dim, tagbase);
          break;
        case OMEGA_H_I32:
          transfer_copy_tmpl<I32>(new_mesh, prod_dim, tagbase);
          break;
        case OMEGA_H_I64:
          transfer_copy_tmpl<I64>(new_mesh, prod_dim, tagbase);
          break;
        case OMEGA_H_F64:
          transfer_copy_tmpl<Real>(new_mesh, prod_dim, tagbase);
          break;
      }
    }
  }
  if (opts.user_xfer) opts.user_xfer->swap_copy_verts(*old_mesh, *new_mesh);
}

template <typename T>
static void transfer_inherit_swap_tmpl(Mesh* old_mesh, Mesh* new_mesh,
    Int prod_dim, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents, TagBase const* tagbase) {
  auto const& name = tagbase->name();
  auto old_tag = old_mesh->get_tag<T>(EDGE, name);
  auto ncomps = old_tag->ncomps();
  auto edge_data = old_tag->array();
  auto key_data = read(unmap(keys2edges, edge_data, ncomps));
  auto prod_data = expand(key_data, keys2prods, ncomps);
  transfer_common(old_mesh, new_mesh, prod_dim, same_ents2old_ents,
      same_ents2new_ents, prods2new_ents, old_tag, prod_data);
}

static void transfer_inherit_swap(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, Int prod_dim, LOs keys2edges, LOs keys2prods,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents) {
  for (Int i = 0; i < old_mesh->ntags(prod_dim); ++i) {
    auto tagbase = old_mesh->get_tag(prod_dim, i);
    if (should_inherit(old_mesh, opts, prod_dim, tagbase)) {
      switch (tagbase->type()) {
        case OMEGA_H_I8:
          transfer_inherit_swap_tmpl<I8>(old_mesh, new_mesh, prod_dim,
              keys2edges, keys2prods, prods2new_ents, same_ents2old_ents,
              same_ents2new_ents, tagbase);
          break;
        case OMEGA_H_I32:
          transfer_inherit_swap_tmpl<I32>(old_mesh, new_mesh, prod_dim,
              keys2edges, keys2prods, prods2new_ents, same_ents2old_ents,
              same_ents2new_ents, tagbase);
          break;
        case OMEGA_H_I64:
          transfer_inherit_swap_tmpl<I64>(old_mesh, new_mesh, prod_dim,
              keys2edges, keys2prods, prods2new_ents, same_ents2old_ents,
              same_ents2new_ents, tagbase);
          break;
        case OMEGA_H_F64:
          transfer_inherit_swap_tmpl<Real>(old_mesh, new_mesh, prod_dim,
              keys2edges, keys2prods, prods2new_ents, same_ents2old_ents,
              same_ents2new_ents, tagbase);
          break;
      }
    }
  }
}

void transfer_swap(Mesh* old_mesh, TransferOpts const& opts, Mesh* new_mesh,
    Int prod_dim, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  begin_code("transfer_swap");
  OMEGA_H_CHECK(prod_dim != VERT);
  transfer_inherit_swap(old_mesh, opts, new_mesh, prod_dim, keys2edges,
      keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents);
  if (prod_dim == EDGE) {
    transfer_length(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
  } else if (prod_dim == FACE) {
    transfer_face_flux(old_mesh, new_mesh, same_ents2old_ents,
        same_ents2new_ents, prods2new_ents);
  }
  if (prod_dim == old_mesh->dim()) {
    transfer_size(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
    transfer_quality(old_mesh, new_mesh, same_ents2old_ents, same_ents2new_ents,
        prods2new_ents);
    transfer_pointwise(old_mesh, opts, new_mesh, EDGE, keys2edges, keys2prods,
        prods2new_ents, same_ents2old_ents, same_ents2new_ents);
    transfer_densities_and_conserve_swap(old_mesh, opts, new_mesh, keys2edges,
        keys2prods, prods2new_ents, same_ents2old_ents, same_ents2new_ents);
  }
  if (opts.user_xfer) {
    opts.user_xfer->swap(*old_mesh, *new_mesh, prod_dim, keys2edges, keys2prods,
        prods2new_ents, same_ents2old_ents, same_ents2new_ents);
  }
  end_code();
}

#define INST(T)                                                                \
  template void transfer_common3(                                              \
      Mesh* new_mesh, Int ent_dim, TagBase const* tagbase, Write<T> new_data); \
  template void transfer_common2(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,  \
      LOs same_ents2old_ents, LOs same_ents2new_ents, TagBase const* tagbase,  \
      Write<T> new_data);                                                      \
  template void transfer_common(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,   \
      LOs same_ents2old_ents, LOs same_ents2new_ents, LOs prods2new_ents,      \
      TagBase const* tagbase, Read<T> prod_data);                              \
  template void transfer_inherit_refine<T>(Mesh * old_mesh, Mesh * new_mesh,   \
      LOs keys2edges, Int prod_dim, LOs keys2prods, LOs prods2new_ents,        \
      LOs same_ents2old_ents, LOs same_ents2new_ents,                          \
      std::string const& name);
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
