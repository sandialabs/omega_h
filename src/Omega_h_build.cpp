#include "Omega_h_build.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_box.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_owners.hpp"
#include "Omega_h_simplify.hpp"

namespace Omega_h {

void add_ents2verts(
    Mesh* mesh, Int ent_dim, LOs ev2v, GOs vert_globals, GOs elem_globals) {
  auto comm = mesh->comm();
  auto deg = ent_dim + 1;
  auto ne = divide_no_remainder(ev2v.size(), deg);
  Remotes owners;
  if (comm->size() > 1) {
    if (mesh->could_be_shared(ent_dim)) {
      if (ent_dim == mesh->dim()) {
        owners = owners_from_globals(comm, elem_globals, Read<I32>());
      } else {
        resolve_derived_copies(comm, vert_globals, deg, &ev2v, &owners);
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
    auto down = reflect_down(ev2v, lv2v, v2l, ent_dim, ldim);
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
    auto mv2v = find_unique(ev2v, elem_dim, mdim);
    add_ents2verts(mesh, mdim, mv2v, vert_globals, elem_globals);
  }
  add_ents2verts(mesh, elem_dim, ev2v, vert_globals, elem_globals);
  if (!comm->reduce_and(is_sorted(vert_globals))) {
    reorder_by_globals(mesh);
  }
}

void build_from_elems2verts(
    Mesh* mesh, CommPtr comm, Int edim, LOs ev2v, Read<GO> vert_globals) {
  mesh->set_comm(comm);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->set_dim(edim);
  build_verts_from_globals(mesh, vert_globals);
  build_ents_from_elems2verts(mesh, ev2v, vert_globals);
}

void build_from_elems2verts(Mesh* mesh, Int edim, LOs ev2v, LO nverts) {
  auto vert_globals = Read<GO>(nverts, 0, 1);
  build_from_elems2verts(
      mesh, mesh->library()->self(), edim, ev2v, vert_globals);
}

void build_from_elems_and_coords(Mesh* mesh, Int edim, LOs ev2v, Reals coords) {
  auto nverts = coords.size() / edim;
  build_from_elems2verts(mesh, edim, ev2v, nverts);
  mesh->add_coords(coords);
}

void build_box_internal(
    Mesh* mesh, Real x, Real y, Real z, LO nx, LO ny, LO nz) {
  OMEGA_H_CHECK(nx > 0);
  OMEGA_H_CHECK(ny >= 0);
  OMEGA_H_CHECK(nz >= 0);
  if (ny == 0) {
    LOs ev2v;
    Reals coords;
    make_1d_box(x, nx, &ev2v, &coords);
    build_from_elems_and_coords(mesh, EDGE, ev2v, coords);
  } else if (nz == 0) {
    LOs qv2v;
    Reals coords;
    make_2d_box(x, y, nx, ny, &qv2v, &coords);
    auto tv2v = simplify::tris_from_quads(qv2v);
    build_from_elems_and_coords(mesh, TRI, tv2v, coords);
  } else {
    LOs hv2v;
    Reals coords;
    make_3d_box(x, y, z, nx, ny, nz, &hv2v, &coords);
    auto tv2v = simplify::tets_from_hexes(hv2v);
    build_from_elems_and_coords(mesh, TET, tv2v, coords);
  }
}

Mesh build_box(CommPtr comm, Real x, Real y, Real z, LO nx, LO ny, LO nz) {
  auto lib = comm->library();
  auto mesh = Mesh(lib);
  if (comm->rank() == 0) {
    build_box_internal(&mesh, x, y, z, nx, ny, nz);
    reorder_by_hilbert(&mesh);
    classify_by_angles(&mesh, PI / 4);
    set_box_class_ids(&mesh, x, y, z, nx, ny, nz);
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
  auto ev2vg = unmap(ev2v, verts2globs, 1);
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
  find_matches_ex(deg, se2fsv, sev2vg, sev2vg, sv2se, &se2ose, &se2ose_codes);
  auto ose2oe = out_dist.items2dests();
  auto se2oe = unmap(se2ose, ose2oe);
  out_dist.set_roots2items(LOs());
  auto e2oe = out_dist.exch(se2oe, 1);
  *p_ents2owners = e2oe;
}

}  // end namespace Omega_h
