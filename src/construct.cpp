#include "construct.hpp"

#include "adjacency.hpp"
#include "align.hpp"
#include "array.hpp"
#include "box.hpp"
#include "linpart.hpp"
#include "map.hpp"
#include "owners.hpp"
#include "remotes.hpp"
#include "simplify.hpp"

namespace osh {

void add_ents2verts(Mesh* mesh, Int edim, LOs ev2v, Read<GO> vert_globals) {
  auto comm = mesh->comm();
  Remotes owners;
  if (comm->size() > 1) {
    auto deg = edim + 1;
    if (edim < mesh->dim()) {
      resolve_derived_copies(comm, vert_globals, deg, &ev2v, &owners);
    } else {
      auto ne = ev2v.size() / deg;
      owners = identity_remotes(comm, ne);
    }
  }
  if (edim == 1) {
    mesh->set_ents(edim, Adj(ev2v));
  } else {
    auto ldim = edim - 1;
    auto lv2v = mesh->ask_verts_of(ldim);
    auto v2l = mesh->ask_up(VERT, ldim);
    auto down = reflect_down(ev2v, lv2v, v2l, edim, ldim);
    mesh->set_ents(edim, down);
  }
  if (comm->size() > 1) {
    mesh->set_owners(edim, owners);
    globals_from_owners(mesh, edim);
  } else {
    mesh->ask_globals(edim);
  }
}

void build_from_elems2verts(
    Mesh* mesh, CommPtr comm, Int edim, LOs ev2v, Read<GO> vert_globals) {
  mesh->set_comm(comm);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->set_dim(edim);
  auto nverts = vert_globals.size();
  mesh->set_verts(nverts);
  mesh->add_tag(VERT, "global", 1, OMEGA_H_GLOBAL, vert_globals);
  if (comm->size() > 1) {
    mesh->set_owners(
        VERT, owners_from_globals(comm, vert_globals, Read<I32>()));
  }
  for (Int mdim = 1; mdim < edim; ++mdim) {
    auto mv2v = find_unique(ev2v, edim, mdim);
    add_ents2verts(mesh, mdim, mv2v, vert_globals);
  }
  add_ents2verts(mesh, edim, ev2v, vert_globals);
}

void build_from_elems2verts(
    Mesh* mesh, Library const& lib, Int edim, LOs ev2v, LO nverts) {
  build_from_elems2verts(mesh, lib.self(), edim, ev2v, Read<GO>(nverts, 0, 1));
}

void build_from_elems_and_coords(
    Mesh* mesh, Library const& lib, Int edim, LOs ev2v, Reals coords) {
  auto nverts = coords.size() / edim;
  build_from_elems2verts(mesh, lib, edim, ev2v, nverts);
  mesh->add_coords(coords);
}

void build_box(Mesh* mesh, Library const& lib, Real x, Real y, Real z, LO nx,
    LO ny, LO nz) {
  CHECK(nx > 0);
  CHECK(ny > 0);
  CHECK(nz >= 0);
  if (nz == 0) {
    LOs qv2v;
    Reals coords;
    make_2d_box(x, y, nx, ny, &qv2v, &coords);
    auto tv2v = simplify::tris_from_quads(qv2v);
    build_from_elems_and_coords(mesh, lib, TRI, tv2v, coords);
  } else {
    LOs hv2v;
    Reals coords;
    make_3d_box(x, y, z, nx, ny, nz, &hv2v, &coords);
    auto tv2v = simplify::tets_from_hexes(hv2v);
    build_from_elems_and_coords(mesh, lib, TET, tv2v, coords);
  }
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

}  // end namespace osh
