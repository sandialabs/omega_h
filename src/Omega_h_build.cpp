#include "Omega_h_build.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_box.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_inertia.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_migrate.hpp"
#include "Omega_h_owners.hpp"
#include "Omega_h_simplify.hpp"

namespace Omega_h {

void add_ents2verts(
    Mesh* mesh, Int ent_dim, LOs ev2v, GOs vert_globals, GOs elem_globals) {
  auto comm = mesh->comm();
  auto nverts_per_ent = element_degree(mesh->family(), ent_dim, VERT);
  auto ne = divide_no_remainder(ev2v.size(), nverts_per_ent);
  Remotes owners;
  if (comm->size() > 1) {
    if (mesh->could_be_shared(ent_dim)) {
      if (ent_dim == mesh->dim()) {
        owners = owners_from_globals(comm, elem_globals, Read<I32>());
      } else {
        resolve_derived_copies(
            comm, vert_globals, nverts_per_ent, &ev2v, &owners);
      }
    } else {
      owners = identity_remotes(comm, ne);
    }
  }
  if (ent_dim == 1) {
    mesh->set_ents(ent_dim, Adj(ev2v));
  } else {
    auto ldim = ent_dim - 1;
    auto lv2v = mesh->ask_verts_of(ldim);
    auto v2l = mesh->ask_up(VERT, ldim);
    auto down = reflect_down(ev2v, lv2v, v2l, mesh->family(), ent_dim, ldim);
    mesh->set_ents(ent_dim, down);
  }
  if (comm->size() > 1) {
    mesh->set_owners(ent_dim, owners);
    if (ent_dim == mesh->dim() && elem_globals.exists()) {
      mesh->add_tag(ent_dim, "global", 1, elem_globals);
    } else {
      globals_from_owners(mesh, ent_dim);
    }
  } else {
    mesh->add_tag(ent_dim, "global", 1, GOs(ne, 0, 1));
  }
}

void build_verts_from_globals(Mesh* mesh, GOs vert_globals) {
  auto comm = mesh->comm();
  auto nverts = vert_globals.size();
  mesh->set_verts(nverts);
  mesh->add_tag(VERT, "global", 1, vert_globals);
  if (comm->size() > 1) {
    mesh->set_owners(
        VERT, owners_from_globals(comm, vert_globals, Read<I32>()));
  }
}

void build_ents_from_elems2verts(
    Mesh* mesh, LOs ev2v, GOs vert_globals, GOs elem_globals) {
  auto comm = mesh->comm();
  auto elem_dim = mesh->dim();
  for (Int mdim = 1; mdim < elem_dim; ++mdim) {
    auto mv2v = find_unique(ev2v, mesh->family(), elem_dim, mdim);
    add_ents2verts(mesh, mdim, mv2v, vert_globals, elem_globals);
  }
  add_ents2verts(mesh, elem_dim, ev2v, vert_globals, elem_globals);
  if (!comm->reduce_and(is_sorted(vert_globals))) {
    reorder_by_globals(mesh);
  }
}

void build_from_elems2verts(Mesh* mesh, CommPtr comm, Omega_h_Family family,
    Int edim, LOs ev2v, Read<GO> vert_globals) {
  mesh->set_comm(comm);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->set_family(family);
  mesh->set_dim(edim);
  build_verts_from_globals(mesh, vert_globals);
  build_ents_from_elems2verts(mesh, ev2v, vert_globals);
}

void build_from_elems2verts(
    Mesh* mesh, Omega_h_Family family, Int edim, LOs ev2v, LO nverts) {
  auto vert_globals = Read<GO>(nverts, 0, 1);
  build_from_elems2verts(
      mesh, mesh->library()->self(), family, edim, ev2v, vert_globals);
}

void build_from_elems_and_coords(
    Mesh* mesh, Omega_h_Family family, Int edim, LOs ev2v, Reals coords) {
  auto nverts = coords.size() / edim;
  build_from_elems2verts(mesh, family, edim, ev2v, nverts);
  mesh->add_coords(coords);
}

void build_box_internal(Mesh* mesh, Omega_h_Family family, Real x, Real y,
    Real z, LO nx, LO ny, LO nz, bool symmetric) {
  OMEGA_H_CHECK(nx > 0);
  OMEGA_H_CHECK(ny >= 0);
  OMEGA_H_CHECK(nz >= 0);
  if (ny == 0) {
    LOs ev2v;
    Reals coords;
    make_1d_box(x, nx, &ev2v, &coords);
    build_from_elems_and_coords(mesh, family, EDGE, ev2v, coords);
  } else if (nz == 0) {
    LOs fv2v;
    Reals coords;
    make_2d_box(x, y, nx, ny, &fv2v, &coords);
    if (family == OMEGA_H_SIMPLEX && (!symmetric)) fv2v = tris_from_quads(fv2v);
    auto fam2 = symmetric ? OMEGA_H_HYPERCUBE : family;
    build_from_elems_and_coords(mesh, fam2, FACE, fv2v, coords);
    if (family == OMEGA_H_SIMPLEX && symmetric) tris_from_quads_symmetric(mesh);
  } else {
    LOs rv2v;
    Reals coords;
    make_3d_box(x, y, z, nx, ny, nz, &rv2v, &coords);
    if (family == OMEGA_H_SIMPLEX && (!symmetric)) rv2v = tets_from_hexes(rv2v);
    auto fam2 = symmetric ? OMEGA_H_HYPERCUBE : family;
    build_from_elems_and_coords(mesh, fam2, REGION, rv2v, coords);
    if (family == OMEGA_H_SIMPLEX && symmetric) tets_from_hexes_symmetric(mesh);
  }
}

Mesh build_box(CommPtr comm, Omega_h_Family family, Real x, Real y, Real z,
    LO nx, LO ny, LO nz, bool symmetric) {
  auto lib = comm->library();
  auto mesh = Mesh(lib);
  if (comm->rank() == 0) {
    build_box_internal(&mesh, family, x, y, z, nx, ny, nz, symmetric);
    reorder_by_hilbert(&mesh);
    classify_box(&mesh, x, y, z, nx, ny, nz);
    mesh.class_sets = get_box_class_sets(mesh.dim());
  }
  mesh.set_comm(comm);
  mesh.balance();
  return mesh;
}

/* When we try to build a mesh from _partitioned_
   element-to-vertex connectivity only, we have to derive
   consistent edges and faces in parallel.
   We'll start by using the usual local derivation on each
   part, which creates the right entities but with inconsistent
   alignment and no ownership information.
   This function establishes the ownership of derived entities
   and canonicalizes their connectivity.
   It does this by expressing connectivity in terms of vertex
   global numbers, locally sorting, and sending each entity
   to its lowest-global-number vertex.
   It uses a linear partitioning of the vertices by global number,
   and each "server" vertex handles the entity copies for which
   it is the lowest-global-number vertex.
   We reuse the matching code that powers reflect_down() to identify
   copies which have the exact same connectivity and establish ownership.
*/

LOs sort_locally_based_on_rank(
    LOs servers_to_served, Read<I32> served_to_rank) {
  OMEGA_H_TIME_FUNCTION;
  auto const served_order = Write<LO>(served_to_rank.size());
  auto functor = OMEGA_H_LAMBDA(LO const server) {
    auto const begin = servers_to_served[server];
    auto const end = servers_to_served[server + 1];
    I32 last_smallest_rank = -1;
    for (LO i = begin; i < end;) {
      I32 next_smallest_rank = ArithTraits<I32>::max();
      for (LO j = begin; j < end; ++j) {
        auto const rank = served_to_rank[j];
        if (rank > last_smallest_rank && rank < next_smallest_rank) {
          next_smallest_rank = rank;
        }
      }
      OMEGA_H_CHECK(next_smallest_rank > last_smallest_rank);
      OMEGA_H_CHECK(next_smallest_rank < ArithTraits<I32>::max());
      for (LO j = begin; j < end; ++j) {
        if (served_to_rank[j] == next_smallest_rank) {
          served_order[i++] = j;
        }
      }
      last_smallest_rank = next_smallest_rank;
    }
  };
  parallel_for(servers_to_served.size() - 1, std::move(functor));
  return served_order;
}

void resolve_derived_copies(CommPtr comm, Read<GO> verts2globs, Int deg,
    LOs* p_ent_verts2verts, Remotes* p_ents2owners) {
  // entity vertices to vertices
  auto const ev2v = *p_ent_verts2verts;
  // entity vertices to vertex globals
  auto const ev2vg = read(unmap(ev2v, verts2globs, 1));
  auto const canon_codes = get_codes_to_canonical(deg, ev2vg);
  // entity vertices to vertices, flipped so smallest global goes in front
  auto const ev2v_canon = align_ev2v(deg, ev2v, canon_codes);
  *p_ent_verts2verts = ev2v_canon;
  // entity vertices to vertex globals, with smallest global first
  auto const ev2vg_canon = align_ev2v(deg, ev2vg, canon_codes);
  // entity to adj vertex with smallest global
  auto const e2fv = get_component(ev2v_canon, deg, 0);
  auto const total_verts = find_total_globals(comm, verts2globs);
  // vertices to owners in linear partitioning of globals
  auto const v2ov = globals_to_linear_owners(comm, verts2globs, total_verts);
  // entities to owning vertices in linear partitioning of vertex globals
  auto const e2ov = unmap(e2fv, v2ov);
  auto const linsize = linear_partition_size(comm, total_verts);
  // dist from entities to linear partitioning of vertex globals
  auto const in_dist = Dist(comm, e2ov, linsize);
  // each MPI rank in the linear partitioning of vertex globals
  // is "serving" a contiguous subset of vertex globals, and hence
  // a subset of the vertices
  // Each of these vertices is "serving" all entities for which they
  // are the smallest-global adjacent vertex
  //
  // "served" entity vertices to vertex globals
  auto const sev2vg = in_dist.exch(ev2vg_canon, deg);
  auto out_dist = in_dist.invert();
  // "serving" vertices to served entities
  auto const sv2svse = out_dist.roots2items();
  // number of served entities
  auto const nse = out_dist.nitems();
  // items2dests() returns, for each served entity, the (rank, local index) pair
  // of where it originally came from
  //
  // owning served entity to original entity
  auto const se2orig = out_dist.items2dests();
  auto const se2orig_rank = se2orig.ranks;
  // the ordering of this array is critical for the resulting output to be
  // both deterministic and serial-parallel consistent.
  // it ensures that entities will be explored for matches in an order that
  // depends only on which rank they came from (which is unique), therefore
  // basing ownership on a minimum-rank rule.
  // note that otherwise the ordering of entities in these arrays is not even
  // deterministic, let alone partition-independent
  auto const svse2se = sort_locally_based_on_rank(sv2svse, se2orig_rank);
  // helper to make "serving vertices to served entities" look like an Adj
  auto const sv2se_codes = Read<I8>(nse, make_code(false, 0, 0));
  auto const sv2se = Adj(sv2svse, svse2se, sv2se_codes);
  // served entities to serving vertex
  auto const se2fsv = invert_fan(sv2svse);
  // find_matches_ex will match up all duplicate entities by assigning a
  // temporary owner to each set of duplicates. note that this owner is entirely
  // ordering-dependent, i.e. it will be just the first duplicate in the input
  // arrays, and the ordering of those arrays is not even deterministic!
  //
  // served entity to owning served entity
  Write<LO> se2ose;
  Write<I8> ignored_codes;
  constexpr bool allow_duplicates = true;
  find_matches_ex(deg, se2fsv, sev2vg, sev2vg, sv2se, &se2ose, &ignored_codes,
      allow_duplicates);
  // the problem with the temporary owner is that it is not deterministic, and
  // definitely not serial-parallel consistent. we need to choose an owner based
  // on a deterministic and serial-parallel consistent rule.
  // Two commonly used rules are:
  //   1) smallest rank
  //   2) rank who owns the fewest entities (tiebreaker is rank).
  //      The problem with this rule is that we are inside the code that
  //      determines ownership, so there is a bit of a chicken-and-egg obstacle
  //      to applying this rule.
  //
  //      TODO: fix non-determinism and partition-dependence of owner
  //
  // served entity to owning original entity
  auto const se2owner_orig = unmap(read(se2ose), se2orig);
  // remove roots2items so that the next exch() call will accept items, not
  // roots
  out_dist.set_roots2items(LOs());
  // send back, to each served entity, the (rank, local index) pair of its new
  // owner
  //
  // entities to owning entities
  auto const e2owner = out_dist.exch(se2owner_orig, 1);
  *p_ents2owners = e2owner;
}

void suggest_slices(
    GO total, I32 comm_size, I32 comm_rank, GO* p_begin, GO* p_end) {
  auto comm_size_gt = GO(comm_size);
  auto quot = total / comm_size_gt;
  auto rem = total % comm_size_gt;
  if (comm_rank < rem) {
    *p_begin = quot * comm_rank + comm_rank;
    *p_end = *p_begin + quot + 1;
  } else {
    *p_begin = quot * comm_rank + rem;
    *p_end = *p_begin + quot;
  }
}

void assemble_slices(CommPtr comm, Omega_h_Family family, Int dim,
    GO global_nelems, GO elem_offset, GOs conn_in, GO global_nverts,
    GO vert_offset, Reals vert_coords, Dist* p_slice_elems2elems,
    LOs* p_conn_out, Dist* p_slice_verts2verts) {
  auto verts_per_elem = element_degree(family, dim, 0);
  auto nslice_elems = divide_no_remainder(conn_in.size(), verts_per_elem);
  auto nslice_verts = divide_no_remainder(vert_coords.size(), dim);
  {  // require that slicing was done as suggested
    GO suggested_elem_offset, suggested_elem_end;
    suggest_slices(global_nelems, comm->size(), comm->rank(),
        &suggested_elem_offset, &suggested_elem_end);
    OMEGA_H_CHECK(suggested_elem_offset == elem_offset);
    OMEGA_H_CHECK(suggested_elem_end == elem_offset + nslice_elems);
    GO suggested_vert_offset, suggested_vert_end;
    suggest_slices(global_nverts, comm->size(), comm->rank(),
        &suggested_vert_offset, &suggested_vert_end);
    OMEGA_H_CHECK(suggested_vert_offset == vert_offset);
    OMEGA_H_CHECK(suggested_vert_end == vert_offset + nslice_verts);
  }
  // generate communication pattern from sliced vertices to their sliced
  // elements
  auto slice_elems2slice_verts = copies_to_linear_owners(comm, conn_in);
  auto slice_verts2slice_elems = slice_elems2slice_verts.invert();
  // compute sliced element coordinates
  auto slice_elem_vert_coords = slice_verts2slice_elems.exch(vert_coords, dim);
  auto slice_elem_coords_w = Write<Real>(nslice_elems * dim);
  auto f = OMEGA_H_LAMBDA(LO slice_elem) {
    for (Int i = 0; i < dim; ++i) {
      slice_elem_coords_w[slice_elem * dim + i] = 0.;
      for (Int j = 0; j < verts_per_elem; ++j) {
        slice_elem_coords_w[slice_elem * dim + i] +=
            slice_elem_vert_coords[(slice_elem * verts_per_elem + j) * dim + i];
      }
      slice_elem_coords_w[slice_elem * dim + i] /= verts_per_elem;
    }
  };
  parallel_for(nslice_elems, f, "sliced elem coords");
  auto slice_elem_coords = Reals(slice_elem_coords_w);
  // geometrically partition the elements (RIB)
  auto rib_masses = Reals(nslice_elems, 1.0);
  auto rib_tol = 2.0;
  auto elem_owners = identity_remotes(comm, nslice_elems);
  auto rib_coords = resize_vectors(slice_elem_coords, dim, 3);
  inertia::Rib rib_hints;
  inertia::recursively_bisect(
      comm, rib_tol, &rib_coords, &rib_masses, &elem_owners, &rib_hints);
  // communication pattern from sliced elements to partitioned elements
  auto elems2slice_elems = Dist{comm, elem_owners, nslice_elems};
  auto slice_elems2elems = elems2slice_elems.invert();
  // communication pattern from sliced vertices to partitioned elements
  auto slice_elem_vert_owners = slice_elems2slice_verts.items2dests();
  auto elem_vert_owners =
      slice_elems2elems.exch(slice_elem_vert_owners, verts_per_elem);
  auto elems2slice_verts = Dist{comm, elem_vert_owners, nslice_verts};
  // unique set of vertices needed for partitioned elements
  auto slice_vert_globals =
      GOs{nslice_verts, vert_offset, 1, "slice vert globals"};
  auto verts2slice_verts =
      get_new_copies2old_owners(elems2slice_verts, slice_vert_globals);
  // new (local) connectivity
  auto slice_verts2elems = elems2slice_verts.invert();
  auto new_conn = form_new_conn(verts2slice_verts, slice_verts2elems);
  *p_slice_elems2elems = slice_elems2elems;
  *p_conn_out = new_conn;
  *p_slice_verts2verts = verts2slice_verts.invert();
}

}  // end namespace Omega_h
