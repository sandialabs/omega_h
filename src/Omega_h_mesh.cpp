#include "Omega_h_mesh.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_bcast.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_control.hpp"
#include "Omega_h_ghost.hpp"
#include "Omega_h_inertia.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_migrate.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_reorder.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_simplex.hpp"
#include "Omega_h_timer.hpp"

namespace Omega_h {

Mesh::Mesh(Library* library) {
  dim_ = -1;
  for (Int i = 0; i <= 3; ++i) nents_[i] = -1;
  parting_ = -1;
  nghost_layers_ = -1;
  OMEGA_H_CHECK(library != nullptr);
  library_ = library;
}

Library* Mesh::library() const { return library_; }

void Mesh::set_comm(CommPtr const& new_comm) {
  auto rank_had_comm = bool(comm_);
  auto nnew_had_comm = new_comm->allreduce(I32(rank_had_comm), OMEGA_H_SUM);
  if (0 < nnew_had_comm && nnew_had_comm < new_comm->size()) {
    // partitioning out from small sub-communicator to larger one
    if (!rank_had_comm) {
      // temporarily set the uninit ranks to Comm::self()
      comm_ = library_->self();
    } else {
      /* forget RIB hints. this prevents some ranks from
         having hints while the new ranks do not, which would
         break RIB. also, since repartitioning does not change
         the geometric properties of the mesh and our RIB code
         is partition and order independent, it will recover
         the same axes of separation as before */
      rib_hints_ = RibPtr();
    }
    bcast_mesh(this, new_comm, rank_had_comm);
  }
  /* if some ranks already have mesh data, their
     parallel info needs updating, we'll do this
     by using the old Dist to set new owners */
  if (0 < nnew_had_comm && library_->world()->size() > 1) {
    for (Int d = 0; d <= dim(); ++d) {
      auto dist = ask_dist(d);
      dist.change_comm(new_comm);
      owners_[d].ranks = dist.items2ranks();
    }
  }
  comm_ = new_comm;
}

void Mesh::set_dim(Int dim) {
  OMEGA_H_CHECK(dim_ == -1);
  OMEGA_H_CHECK(dim >= 1);
  OMEGA_H_CHECK(dim <= 3);
  dim_ = dim;
}

void Mesh::set_verts(LO nverts) { nents_[VERT] = nverts; }

void Mesh::set_ents(Int dim, Adj down) {
  check_dim(dim);
  OMEGA_H_CHECK(!has_ents(dim));
  LOs hl2l = down.ab2b;
  OMEGA_H_CHECK(hl2l.size() % simplex_degrees[dim][dim - 1] == 0);
  nents_[dim] = hl2l.size() / simplex_degrees[dim][dim - 1];
  add_adj(dim, dim - 1, down);
}

CommPtr Mesh::comm() const { return comm_; }

LO Mesh::nents(Int dim) const {
  check_dim2(dim);
  return nents_[dim];
}

LO Mesh::nelems() const { return nents(dim()); }

LO Mesh::ntets() const { return nents(TET); }

LO Mesh::ntris() const { return nents(TRI); }

LO Mesh::nedges() const { return nents(EDGE); }

LO Mesh::nverts() const { return nents(VERT); }

GO Mesh::nglobal_ents(Int dim) {
  if (!could_be_shared(dim)) {
    return comm_->allreduce(GO(nents(dim)), OMEGA_H_SUM);
  }
  auto nowned = get_sum(this->owned(dim));
  return comm_->allreduce(GO(nowned), OMEGA_H_SUM);
}

template <typename T>
void Mesh::add_tag(Int dim, std::string const& name, Int ncomps) {
  check_dim2(dim);
  check_tag_name(name);
  if (has_tag(dim, name)) {
    Omega_h_fail(
        "add_tag(%s, %s): already exists. use set_tag or "
        "remove_tag\n",
        plural_names[dim], name.c_str());
  }
  OMEGA_H_CHECK(ncomps >= 0);
  OMEGA_H_CHECK(ncomps <= Int(INT8_MAX));
  OMEGA_H_CHECK(tags_[dim].size() < size_t(INT8_MAX));
  tags_[dim].push_back(TagPtr(new Tag<T>(name, ncomps)));
}

template <typename T>
void Mesh::add_tag(Int dim, std::string const& name, Int ncomps, Read<T> array,
    bool internal) {
  add_tag<T>(dim, name, ncomps);
  set_tag<T>(dim, name, array, internal);
}

template <typename T>
void Mesh::set_tag(
    Int dim, std::string const& name, Read<T> array, bool internal) {
  if (!has_tag(dim, name)) {
    Omega_h_fail("set_tag(%s, %s): tag doesn't exist (use add_tag first)\n",
        plural_names[dim], name.c_str());
  }
  Tag<T>* tag = as<T>(tag_iter(dim, name)->get());
  OMEGA_H_CHECK(array.size() == nents(dim) * tag->ncomps());
  /* internal typically indicates migration/adaptation/file reading,
     when we do not want any invalidation to take place.
     the invalidation is there to prevent users changing coordinates
     etc. without updating dependent fields */
  if (!internal) react_to_set_tag(dim, name);
  tag->set_array(array);
}

void Mesh::react_to_set_tag(Int dim, std::string const& name) {
  /* hardcoded cache invalidations */
  if ((dim == VERT) && ((name == "coordinates") || (name == "metric"))) {
    remove_tag(EDGE, "length");
    remove_tag(this->dim(), "quality");
  }
  if ((dim == VERT) && (name == "coordinates")) {
    remove_tag(this->dim(), "size");
  }
}

TagBase const* Mesh::get_tagbase(Int dim, std::string const& name) const {
  check_dim2(dim);
  if (!has_tag(dim, name)) {
    Omega_h_fail("get_tagbase(%s, %s): doesn't exist\n", plural_names[dim],
        name.c_str());
  }
  return tag_iter(dim, name)->get();
}

template <typename T>
Tag<T> const* Mesh::get_tag(Int dim, std::string const& name) const {
  return as<T>(get_tagbase(dim, name));
}

template <typename T>
Read<T> Mesh::get_array(Int dim, std::string const& name) const {
  return get_tag<T>(dim, name)->array();
}

void Mesh::remove_tag(Int dim, std::string const& name) {
  if (!has_tag(dim, name)) return;
  check_dim2(dim);
  OMEGA_H_CHECK(has_tag(dim, name));
  tags_[dim].erase(tag_iter(dim, name));
}

bool Mesh::has_tag(Int dim, std::string const& name) const {
  check_dim(dim);
  if (!has_ents(dim)) return false;
  return tag_iter(dim, name) != tags_[dim].end();
}

Int Mesh::ntags(Int dim) const {
  check_dim2(dim);
  return static_cast<Int>(tags_[dim].size());
}

TagBase const* Mesh::get_tag(Int dim, Int i) const {
  check_dim2(dim);
  OMEGA_H_CHECK(0 <= i);
  OMEGA_H_CHECK(i <= ntags(dim));
  return tags_[dim][static_cast<std::size_t>(i)].get();
}

bool Mesh::has_ents(Int dim) const {
  check_dim(dim);
  return nents_[dim] >= 0;
}

bool Mesh::has_adj(Int from, Int to) const {
  check_dim(from);
  check_dim(to);
  return bool(adjs_[from][to]);
}

Adj Mesh::get_adj(Int from, Int to) const {
  check_dim2(from);
  check_dim2(to);
  OMEGA_H_CHECK(has_adj(from, to));
  return *(adjs_[from][to]);
}

Adj Mesh::ask_down(Int from, Int to) {
  OMEGA_H_CHECK(to < from);
  return ask_adj(from, to);
}

LOs Mesh::ask_verts_of(Int dim) { return ask_adj(dim, VERT).ab2b; }

LOs Mesh::ask_elem_verts() { return ask_verts_of(dim()); }

Adj Mesh::ask_up(Int from, Int to) {
  OMEGA_H_CHECK(from < to);
  return ask_adj(from, to);
}

Graph Mesh::ask_star(Int dim) {
  OMEGA_H_CHECK(dim < this->dim());
  return ask_adj(dim, dim);
}

Graph Mesh::ask_dual() { return ask_adj(dim(), dim()); }

struct HasName {
  std::string const& name_;
  HasName(std::string const& name) : name_(name) {}
  bool operator()(std::shared_ptr<TagBase> const& a) {
    return a->name() == name_;
  }
};

Mesh::TagIter Mesh::tag_iter(Int dim, std::string const& name) {
  return std::find_if(tags_[dim].begin(), tags_[dim].end(), HasName(name));
}

Mesh::TagCIter Mesh::tag_iter(Int dim, std::string const& name) const {
  return std::find_if(tags_[dim].begin(), tags_[dim].end(), HasName(name));
}

void Mesh::check_dim(Int dim) const {
  OMEGA_H_CHECK(0 <= dim);
  OMEGA_H_CHECK(dim <= this->dim());
}

void Mesh::check_dim2(Int dim) const {
  check_dim(dim);
  OMEGA_H_CHECK(has_ents(dim));
}

void Mesh::add_adj(Int from, Int to, Adj adj) {
  check_dim2(from);
  check_dim(to);
  OMEGA_H_CHECK(adj.ab2b.exists());
  if (to < from) {
    OMEGA_H_CHECK(!adj.a2ab.exists());
    if (to == VERT) {
      OMEGA_H_CHECK(!adj.codes.exists());
    } else {
      OMEGA_H_CHECK(adj.codes.exists());
    }
    OMEGA_H_CHECK(adj.ab2b.size() == nents(from) * simplex_degrees[from][to]);
  } else {
    if (from < to) {
      OMEGA_H_CHECK(adj.a2ab.exists());
      OMEGA_H_CHECK(adj.codes.exists());
      OMEGA_H_CHECK(adj.ab2b.size() == nents(to) * simplex_degrees[to][from]);
    }
    OMEGA_H_CHECK(adj.a2ab.size() == nents(from) + 1);
  }
  adjs_[from][to] = std::make_shared<Adj>(adj);
}

Adj Mesh::derive_adj(Int from, Int to) {
  check_dim(from);
  check_dim2(to);
  if (from < to) {
    Adj down = ask_adj(to, from);
    Int nlows_per_high = simplex_degrees[to][from];
    LO nlows = nents(from);
    Adj up = invert_adj(down, nlows_per_high, nlows);
    return up;
  } else if (to < from) {
    OMEGA_H_CHECK(to + 1 < from);
    Adj h2m = ask_adj(from, to + 1);
    Adj m2l = ask_adj(to + 1, to);
    Adj h2l = transit(h2m, m2l, from, to);
    return h2l;
  } else {
    if (from == dim() && to == dim()) {
      return elements_across_sides(dim(), ask_adj(dim(), dim() - 1),
          ask_adj(dim() - 1, dim()), mark_exposed_sides(this));
    }
    if (from == VERT && to == VERT) {
      return verts_across_edges(ask_adj(EDGE, VERT), ask_adj(VERT, EDGE));
    }
    if (from == EDGE && to == EDGE) {
      OMEGA_H_CHECK(dim() >= 2);
      Graph g = edges_across_tris(ask_adj(TRI, EDGE), ask_adj(EDGE, TRI));
      if (dim() == 3) {
        g = add_edges(
            g, edges_across_tets(ask_adj(TET, EDGE), ask_adj(EDGE, TET)));
      }
      return g;
    }
  }
  Omega_h_fail("can't derive adjacency from %s to %s\n", plural_names[from],
      plural_names[to]);
  OMEGA_H_NORETURN(Adj());
}

Adj Mesh::ask_adj(Int from, Int to) {
  check_dim2(from);
  check_dim2(to);
  if (has_adj(from, to)) {
    return get_adj(from, to);
  }
  auto t0 = now();
  Adj derived = derive_adj(from, to);
  auto t1 = now();
  add_to_global_timer("deriving adjacencies", t1 - t0);
  adjs_[from][to] = std::make_shared<Adj>(derived);
  return derived;
}

void Mesh::add_coords(Reals array) {
  add_tag<Real>(0, "coordinates", dim(), array);
}

Reals Mesh::coords() const { return get_array<Real>(0, "coordinates"); }

void Mesh::set_coords(Reals const& array) {
  OMEGA_H_CHECK(array.size() == nverts() * dim());
  set_tag<Real>(VERT, "coordinates", array);
}

Read<GO> Mesh::globals(Int dim) const {
  return get_array<GO>(dim, "global");
}

Reals Mesh::ask_lengths() {
  if (!has_tag(EDGE, "length")) {
    auto lengths = measure_edges_metric(this);
    add_tag(EDGE, "length", 1, lengths);
  }
  return get_array<Real>(EDGE, "length");
}

Reals Mesh::ask_qualities() {
  if (!has_tag(dim(), "quality")) {
    auto qualities = measure_qualities(this);
    add_tag(dim(), "quality", 1, qualities);
  }
  return get_array<Real>(dim(), "quality");
}

Reals Mesh::ask_sizes() {
  if (!has_tag(dim(), "size")) {
    auto sizes = measure_elements_real(this);
    add_tag(dim(), "size", 1, sizes);
  }
  return get_array<Real>(dim(), "size");
}

void Mesh::set_owners(Int dim, Remotes owners) {
  check_dim2(dim);
  OMEGA_H_CHECK(nents(dim) == owners.ranks.size());
  OMEGA_H_CHECK(nents(dim) == owners.idxs.size());
  owners_[dim] = owners;
  dists_[dim] = DistPtr();
}

Remotes Mesh::ask_owners(Int dim) {
  if (!owners_[dim].ranks.exists() || !owners_[dim].idxs.exists()) {
    OMEGA_H_CHECK(comm_->size() == 1);
    owners_[dim] =
        Remotes(Read<I32>(nents(dim), comm_->rank()), LOs(nents(dim), 0, 1));
  }
  return owners_[dim];
}

Read<I8> Mesh::owned(Int dim) {
  auto e2rank = ask_owners(dim).ranks;
  return each_eq_to(e2rank, comm()->rank());
}

Dist Mesh::ask_dist(Int dim) {
  if (!dists_[dim]) {
    auto owners = ask_owners(dim);
    OMEGA_H_CHECK(owners.ranks.exists());
    OMEGA_H_CHECK(owners.idxs.exists());
    dists_[dim] = std::make_shared<Dist>(comm_, owners, nents(dim));
  }
  return *(dists_[dim]);
}

Omega_h_Parting Mesh::parting() const {
  OMEGA_H_CHECK(parting_ != -1);
  return Omega_h_Parting(parting_);
}

Int Mesh::nghost_layers() const { return nghost_layers_; }

void Mesh::set_parting(Omega_h_Parting parting, Int nlayers, bool verbose) {
  if (verbose && comm_->rank() == 0) {
    std::cout << "going to ";
    switch (parting) {
      case OMEGA_H_ELEM_BASED:
        std::cout << "element based";
        break;
      case OMEGA_H_VERT_BASED:
        std::cout << "vertex based";
        break;
      case OMEGA_H_GHOSTED:
        std::cout << "ghosted (" << nlayers << " layers)";
        break;
    };
    std::cout << " partitioning\n";
  }
  if (parting_ == -1) {
    parting_ = parting;
    nghost_layers_ = nlayers;
    return;
  }
  if (parting_ == parting && nghost_layers_ == nlayers) {
    return;
  }
  if (parting == OMEGA_H_ELEM_BASED) {
    OMEGA_H_CHECK(nlayers == 0);
    if (comm_->size() > 1) partition_by_elems(this, verbose);
  } else if (parting == OMEGA_H_GHOSTED) {
    if (parting_ != OMEGA_H_GHOSTED || nlayers < nghost_layers_) {
      set_parting(OMEGA_H_ELEM_BASED, 0, false);
    }
    if (comm_->size() > 1) ghost_mesh(this, nlayers, verbose);
  } else if (parting == OMEGA_H_VERT_BASED) {
    OMEGA_H_CHECK(nlayers == 1);
    if (comm_->size() > 1) partition_by_verts(this, verbose);
  }
  parting_ = parting;
  nghost_layers_ = nlayers;
}

void Mesh::set_parting(Omega_h_Parting parting, bool verbose) {
  if (parting == OMEGA_H_ELEM_BASED)
    set_parting(parting, 0, verbose);
  else
    set_parting(parting, 1, verbose);
}

void Mesh::migrate(Remotes new_elems2old_owners, bool verbose) {
  migrate_mesh(this, new_elems2old_owners, verbose);
}

void Mesh::reorder() { reorder_by_hilbert(this); }

void Mesh::balance(bool predictive) {
  if (comm_->size() == 1) return;
  set_parting(OMEGA_H_ELEM_BASED);
  inertia::Rib hints;
  if (rib_hints_) hints = *rib_hints_;
  auto ecoords =
      average_field(this, dim(), LOs(nelems(), 0, 1), dim(), coords());
  if (dim() < 3) ecoords = resize_vectors(ecoords, dim(), 3);
  Reals masses;
  Real abs_tol;
  if (predictive) {
    masses =
        get_expected_nelems_per_elem(this, get_array<Real>(VERT, "metric"));
    /* average between input mesh weight (1.0)
       and predicted output mesh weight */
    masses = add_to_each(masses, 1.);
    masses = multiply_each_by(1. / 2., masses);
    abs_tol = max2(0.0, get_max(comm_, masses));
  } else {
    masses = Reals(nelems(), 1);
    abs_tol = 1.0;
  }
  abs_tol *= 2.0;  // fudge factor ?
  auto owners = ask_owners(dim());
  recursively_bisect(comm(), abs_tol, &ecoords, &masses, &owners, &hints);
  rib_hints_ = std::make_shared<inertia::Rib>(hints);
  migrate(owners);
}

Graph Mesh::ask_graph(Int from, Int to) {
  if (to > from) {
    return ask_up(from, to);
  }
  if (to < from) {
    auto down = ask_down(from, to);
    auto a2ab = LOs(nents(from) + 1, 0, simplex_degrees[from][to]);
    return Graph(a2ab, down.ab2b);
  }
  OMEGA_H_CHECK(from == to);
  return identity_graph(nents(from));
}

template <typename T>
Read<T> Mesh::sync_array(Int ent_dim, Read<T> a, Int width) {
  if (!could_be_shared(ent_dim)) return a;
  return ask_dist(ent_dim).invert().exch(a, width);
}

template <typename T>
Read<T> Mesh::sync_subset_array(
    Int ent_dim, Read<T> a_data, LOs a2e, T default_val, Int width) {
  if (!could_be_shared(ent_dim)) return a_data;
  auto e_data = map_onto(a_data, a2e, nents(ent_dim), default_val, width);
  e_data = sync_array(ent_dim, e_data, width);
  return unmap(a2e, e_data, width);
}

template <typename T>
Read<T> Mesh::reduce_array(Int ent_dim, Read<T> a, Int width, Omega_h_Op op) {
  if (!could_be_shared(ent_dim)) return a;
  return ask_dist(ent_dim).exch_reduce(a, width, op);
}

template <typename T>
Read<T> Mesh::owned_array(Int ent_dim, Read<T> a, Int width) {
  OMEGA_H_CHECK(a.size() == width * nents(ent_dim));
  if (!could_be_shared(ent_dim)) return a;
  auto o = owned(ent_dim);
  auto o2e = collect_marked(o);
  return unmap(o2e, a, width);
}

void Mesh::sync_tag(Int dim, std::string const& name) {
  auto tagbase = get_tagbase(dim, name);
  switch (tagbase->type()) {
    case OMEGA_H_I8: {
      auto out = sync_array(dim, as<I8>(tagbase)->array(), tagbase->ncomps());
      set_tag(dim, name, out);
      break;
    }
    case OMEGA_H_I32: {
      auto out = sync_array(dim, as<I32>(tagbase)->array(), tagbase->ncomps());
      set_tag(dim, name, out);
      break;
    }
    case OMEGA_H_I64: {
      auto out = sync_array(dim, as<I64>(tagbase)->array(), tagbase->ncomps());
      set_tag(dim, name, out);
      break;
    }
    case OMEGA_H_F64: {
      auto out = sync_array(dim, as<Real>(tagbase)->array(), tagbase->ncomps());
      set_tag(dim, name, out);
      break;
    }
  }
}

void Mesh::reduce_tag(Int dim, std::string const& name, Omega_h_Op op) {
  auto tagbase = get_tagbase(dim, name);
  switch (tagbase->type()) {
    case OMEGA_H_I8: {
      auto out =
          reduce_array(dim, as<I8>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(dim, name, out);
      break;
    }
    case OMEGA_H_I32: {
      auto out =
          reduce_array(dim, as<I32>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(dim, name, out);
      break;
    }
    case OMEGA_H_I64: {
      auto out =
          reduce_array(dim, as<I64>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(dim, name, out);
      break;
    }
    case OMEGA_H_F64: {
      auto out =
          reduce_array(dim, as<Real>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(dim, name, out);
      break;
    }
  }
}

bool Mesh::operator==(Mesh& other) {
  auto opts = MeshCompareOpts::init(this, VarCompareOpts::zero_tolerance());
  return OMEGA_H_SAME == compare_meshes(this, &other, opts, false);
}

Real Mesh::min_quality() { return get_min(comm_, ask_qualities()); }

Real Mesh::max_length() { return get_max(comm_, ask_lengths()); }

bool Mesh::could_be_shared(Int ent_dim) const {
  return !((comm_->size() == 1) ||
           (parting_ == OMEGA_H_ELEM_BASED && ent_dim == dim()));
}

bool Mesh::owners_have_all_upward(Int ent_dim) const {
  return have_all_upward() ||
         (parting_ == OMEGA_H_VERT_BASED && ent_dim == VERT);
}

bool Mesh::have_all_upward() const {
  return (comm_->size() == 1) || (parting_ == OMEGA_H_GHOSTED);
}

Mesh Mesh::copy_meta() const {
  Mesh m(library_);
  m.dim_ = this->dim_;
  m.comm_ = this->comm_;
  m.parting_ = this->parting_;
  m.nghost_layers_ = this->nghost_layers_;
  m.rib_hints_ = this->rib_hints_;
  return m;
}

Mesh::RibPtr Mesh::rib_hints() const { return rib_hints_; }

void Mesh::set_rib_hints(RibPtr hints) { rib_hints_ = hints; }

Real Mesh::imbalance(Int ent_dim) const {
  if (ent_dim == -1) ent_dim = this->dim();
  auto local = Real(this->nents(ent_dim));
  auto s = comm_->allreduce(local, OMEGA_H_SUM);
  if (s == 0.0) return 1.0;
  auto m = comm_->allreduce(local, OMEGA_H_MAX);
  auto n = comm_->size();
  auto a = s / n;
  return m / a;
}

bool can_print(Mesh* mesh) {
  return (!mesh->library()->silent_) && (mesh->comm()->rank() == 0);
}

Real repro_sum_owned(Mesh* mesh, Int dim, Reals a) {
  return repro_sum(mesh->comm(), mesh->owned_array(dim, a, 1));
}

Reals average_field(Mesh* mesh, Int dim, LOs a2e, Int ncomps, Reals v2x) {
  auto ev2v = mesh->ask_verts_of(dim);
  auto degree = simplex_degrees[dim][VERT];
  OMEGA_H_CHECK(v2x.size() % ncomps == 0);
  auto na = a2e.size();
  Write<Real> out(na * ncomps);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto e = a2e[a];
    for (Int j = 0; j < ncomps; ++j) {
      Real comp = 0;
      for (Int k = 0; k < degree; ++k) {
        auto v = ev2v[e * degree + k];
        comp += v2x[v * ncomps + j];
      }
      comp /= degree;
      out[a * ncomps + j] = comp;
    }
  };
  parallel_for(na, f);
  return out;
}

Reals average_field(Mesh* mesh, Int dim, Int ncomps, Reals v2x) {
  auto a2e = LOs(mesh->nents(dim), 0, 1);
  return average_field(mesh, dim, a2e, ncomps, v2x);
}

TagSet get_all_mesh_tags(Mesh* mesh) {
  TagSet out;
  for (Int i = 0; i <= mesh->dim(); ++i) {
    for (Int j = 0; j < mesh->ntags(i); ++j) {
      auto tagbase = mesh->get_tag(i, j);
      out[size_t(i)].insert(tagbase->name());
    }
  }
  return out;
}

void ask_for_mesh_tags(Mesh* mesh, TagSet const& tags) {
  if (tags[EDGE].count("length")) mesh->ask_lengths();
  if (tags[size_t(mesh->dim())].count("quality")) mesh->ask_qualities();
}

#define OMEGA_H_INST(T)                                                        \
  template Tag<T> const* Mesh::get_tag<T>(Int dim, std::string const& name)    \
      const;                                                                   \
  template Read<T> Mesh::get_array<T>(Int dim, std::string const& name) const; \
  template void Mesh::add_tag<T>(                                              \
      Int dim, std::string const& name, Int ncomps);                           \
  template void Mesh::add_tag<T>(Int dim, std::string const& name, Int ncomps, \
      Read<T> array, bool internal);                                           \
  template void Mesh::set_tag(                                                 \
      Int dim, std::string const& name, Read<T> array, bool internal);         \
  template Read<T> Mesh::sync_array(Int ent_dim, Read<T> a, Int width);        \
  template Read<T> Mesh::owned_array(Int ent_dim, Read<T> a, Int width);       \
  template Read<T> Mesh::sync_subset_array(                                    \
      Int ent_dim, Read<T> a_data, LOs a2e, T default_val, Int width);         \
  template Read<T> Mesh::reduce_array(                                         \
      Int ent_dim, Read<T> a, Int width, Omega_h_Op op);
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

}  // end namespace Omega_h
