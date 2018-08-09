#include "Omega_h_build.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_box.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_inertia.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_for.hpp"
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

void resolve_derived_copies(CommPtr comm, Read<GO> verts2globs, Int deg,
    LOs* p_ent_verts2verts, Remotes* p_ents2owners) {
  auto ev2v = *p_ent_verts2verts;
  auto ev2vg = read(unmap(ev2v, verts2globs, 1));
  auto canon_codes = get_codes_to_canonical(deg, ev2vg);
  auto ev2v_canon = align_ev2v(deg, ev2v, canon_codes);
  *p_ent_verts2verts = ev2v_canon;
  auto ev2vg_canon = align_ev2v(deg, ev2vg, canon_codes);
  auto e2fv = get_component(ev2v_canon, deg, 0);
  auto total_verts = find_total_globals(comm, verts2globs);
  auto v2ov = globals_to_linear_owners(comm, verts2globs, total_verts);
  auto e2ov = unmap(e2fv, v2ov);
  auto linsize = linear_partition_size(comm, total_verts);
  auto in_dist = Dist(comm, e2ov, linsize);
  auto sev2vg = in_dist.exch(ev2vg_canon, deg);
  auto out_dist = in_dist.invert();
  auto sv2svse = out_dist.roots2items();
  auto nse = out_dist.nitems();
  auto svse2se = LOs(nse, 0, 1);
  auto sv2se_codes = Read<I8>(nse, make_code(false, 0, 0));
  auto sv2se = Adj(sv2svse, svse2se, sv2se_codes);
  auto se2fsv = invert_fan(sv2svse);
  LOs se2ose;
  Read<I8> se2ose_codes;
  constexpr bool allow_duplicates = true;
  find_matches_ex(deg, se2fsv, sev2vg, sev2vg, sv2se, &se2ose, &se2ose_codes,
      allow_duplicates);
  auto ose2oe = out_dist.items2dests();
  auto se2oe = unmap(se2ose, ose2oe);
  out_dist.set_roots2items(LOs());
  auto e2oe = out_dist.exch(se2oe, 1);
  *p_ents2owners = e2oe;
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
