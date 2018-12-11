#include "Omega_h_simplify.hpp"

#include "Omega_h_build.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_hypercube.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_scalar.hpp"

namespace Omega_h {

/* The algorithms here are described in:

   Dompierre, Julien, et al.
   "How to Subdivide Pyramids, Prisms, and Hexahedra into Tetrahedra."
   Proceedings of the 8th International Meshing Roundtable. 1999.
*/

namespace {

OMEGA_H_INLINE Int find_min(LO const v[], Int n) {
  Int min_i = 0;
  LO min_v = v[0];
  for (Int i = 1; i < n; ++i) {
    if (v[i] < min_v) {
      min_i = i;
      min_v = v[i];
    }
  }
  return min_i;
}

OMEGA_H_INLINE bool lt(Int v_i, Int v_j, Int v_k, Int v_l) {
  return min2(v_i, v_j) < min2(v_k, v_l);
}

OMEGA_H_INLINE void rot(LO v[], Int n) {
  auto const tmp = v[n - 1];
  for (Int i = 0; i < n - 1; ++i) v[i + 1] = v[i];
  v[0] = tmp;
}

OMEGA_H_INLINE void rot_ntimes(LO v[], Int nv, Int ntimes) {
  for (Int i = 0; i < ntimes; ++i) rot(v, nv);
}

OMEGA_H_INLINE void rot_to_first(LO v[], Int nv, Int first) {
  rot_ntimes(v, nv, ((nv - first) % nv));
}

/* quad to tri template */
OMEGA_H_CONSTANT_DATA Int const qtv2qqv[2][3] = {{0, 1, 2}, {2, 3, 0}};

/* below are the hex-to-tet templates for the
   four unique cases identified by Dompierre et al.
   They list multiple templates for each case,
   here we'll be a bit lazy and just pick one each */

/* tets from a hex with no diagonals
   into the back-upper-right corner */
OMEGA_H_CONSTANT_DATA Int const htv2hhv_0[5][4] = {
    {0, 1, 2, 5}, {0, 2, 7, 5}, {0, 2, 3, 7}, {0, 5, 7, 4}, {2, 7, 5, 6}};
/* tets from a hex with 1 diagonal
   into the back-upper-right corner,
   on the right face */
OMEGA_H_CONSTANT_DATA Int const htv2hhv_1[6][4] = {{0, 5, 7, 4}, {0, 1, 7, 5},
    {1, 6, 7, 5}, {0, 7, 2, 3}, {0, 7, 1, 2}, {1, 7, 6, 2}};
/* tets from a hex with 2 diagonals
   into the back-upper-right corner,
   none on the right face */
OMEGA_H_CONSTANT_DATA Int const htv2hhv_2[6][4] = {{0, 4, 5, 6}, {0, 3, 7, 6},
    {0, 7, 4, 6}, {0, 1, 2, 5}, {0, 3, 6, 2}, {0, 6, 5, 2}};
/* tets from a hex with 3 diagonals
   into the back-upper-right corner */
OMEGA_H_CONSTANT_DATA Int const htv2hhv_3[6][4] = {{0, 2, 3, 6}, {0, 3, 7, 6},
    {0, 7, 4, 6}, {0, 5, 6, 4}, {1, 5, 6, 0}, {1, 6, 2, 0}};

OMEGA_H_CONSTANT_DATA Int const hex_flip_pairs[4][2] = {
    {0, 4}, {3, 5}, {1, 7}, {2, 6}};

OMEGA_H_DEVICE void flip_hex(LO hhv2v[]) {
  for (Int i = 0; i < 4; ++i)
    swap2(hhv2v[hex_flip_pairs[i][0]], hhv2v[hex_flip_pairs[i][1]]);
}

/* faces adjacent to the back-upper-right
   corner, starting with the right face
   and curling around the centroidal XYZ axis.
   also, numbered with the corner vertex first */
OMEGA_H_CONSTANT_DATA Int const hex_bur_faces[3][4] = {
    {6, 5, 1, 2}, {6, 2, 3, 7}, {6, 7, 4, 5}};

/* the vertices that rotate amonst one another
   when a hex is rotated around its centroidal
   XYZ axis (the line between the front-lower-left
   corner and the back-upper-right corner). */
OMEGA_H_CONSTANT_DATA Int const hex_bur_ring[6] = {1, 2, 3, 7, 5, 4};

OMEGA_H_DEVICE void hex_bur_rot_ntimes(LO hhv2v[], Int ntimes) {
  LO tmp[6];
  for (Int i = 0; i < 6; ++i) {
    tmp[i] = hhv2v[hex_bur_ring[i]];
  }
  rot_ntimes(tmp, 6, ntimes * 2);
  for (Int i = 0; i < 6; ++i) {
    hhv2v[hex_bur_ring[i]] = tmp[i];
  }
}

/* rotate the hex around its centroidal XYZ axis such
   that face (new_right) becomes the right face.
   (new_right) corresponds to the table hex_bur_faces[] */
OMEGA_H_DEVICE void hex_bur_rot_to_right(LO hhv2v[], Int new_right) {
  hex_bur_rot_ntimes(hhv2v, ((3 - new_right) % 3));
}
}  // namespace

LOs tris_from_quads(LOs qv2v) {
  LO nq = divide_no_remainder(qv2v.size(), 4);
  LO nt = nq * 2;
  Write<LO> tv2v(nt * 3);
  auto f = OMEGA_H_LAMBDA(LO q) {
    LO qv_begin = q * 4;
    LO qqv2v[4];
    for (Int i = 0; i < 4; ++i) qqv2v[i] = qv2v[qv_begin + i];
    /* rotate the quad such that the smallest vertex
       is in the lower left */
    Int min_i = find_min(qqv2v, 4);
    rot_to_first(qqv2v, 4, min_i);
    /* split it into two triangles via a template */
    for (Int i = 0; i < 2; ++i) {
      LO t = q * 2 + i;
      LO tv_begin = t * 3;
      for (Int j = 0; j < 3; ++j) {
        tv2v[tv_begin + j] = qqv2v[qtv2qqv[i][j]];
      }
    }
  };
  parallel_for(nq, f, "tris_from_quads");
  return tv2v;
}

static OMEGA_H_DEVICE void tets_from_hex_1(
    LO h, LOs hv2v, LO hhv2v[], Int diags_into[], Int& ndiags_into) {
  LO hv_begin = h * 8;
  for (Int i = 0; i < 8; ++i) {
    hhv2v[i] = hv2v[hv_begin + i];
  }
  /* rotate around the centroidal Z axis
     such that the minimum vertex is the
     front-lower-left or front-upper-left. */
  Int min_i = find_min(hhv2v, 8);
  rot_to_first(hhv2v + 0, 4, min_i % 4);
  rot_to_first(hhv2v + 4, 4, min_i % 4);
  /* rotate around the centroidal XY line
     as needed to make the minimum vertex
     front-lower-left */
  if (min_i >= 4) flip_hex(hhv2v);
  /* check the diagonals going into the
     back-upper-right corner */
  for (Int i = 0; i < 3; ++i) {
    diags_into[i] = lt(hhv2v[hex_bur_faces[i][0]], hhv2v[hex_bur_faces[i][2]],
        hhv2v[hex_bur_faces[i][1]], hhv2v[hex_bur_faces[i][3]]);
  }
  ndiags_into = diags_into[0] + diags_into[1] + diags_into[2];
}

static OMEGA_H_DEVICE void fill_tets_from_hex(Write<LO> tv2v, LOs h2ht, LO h,
    LO const hhv2v[], Int const case_template[][4], Int nhht) {
  LO t = h2ht[h];
  for (Int i = 0; i < nhht; ++i) {
    for (Int j = 0; j < 4; ++j) {
      tv2v[t * 4 + j] = hhv2v[case_template[i][j]];
    }
    ++t;
  }
}

LOs tets_from_hexes(LOs hv2v) {
  LO nh = divide_no_remainder(hv2v.size(), 8);
  Write<LO> degrees(nh);
  auto count = OMEGA_H_LAMBDA(LO h) {
    LO hhv2v[8];
    Int diags_into[3];
    Int ndiags_into;
    tets_from_hex_1(h, hv2v, hhv2v, diags_into, ndiags_into);
    if (ndiags_into == 0)
      degrees[h] = 5;
    else
      degrees[h] = 6;
  };
  parallel_for(nh, count, "tets_from_hexes(count)");
  auto h2ht = offset_scan(LOs(degrees));
  LO nt = h2ht.last();
  Write<LO> tv2v(nt * 4);
  auto fill = OMEGA_H_LAMBDA(LO h) {
    LO hhv2v[8];
    Int diags_into[3];
    Int ndiags_into;
    tets_from_hex_1(h, hv2v, hhv2v, diags_into, ndiags_into);
    if (ndiags_into == 0) {
      fill_tets_from_hex(tv2v, h2ht, h, hhv2v, htv2hhv_0, 5);
    } else if (ndiags_into == 1) {
      /* one diagonal into the far corner.
         find it, rotate it to the right side, apply template */
      Int diag_face = -1;
      for (Int i = 0; i < 3; ++i)
        if (diags_into[i]) diag_face = i;
      hex_bur_rot_to_right(hhv2v, diag_face);
      fill_tets_from_hex(tv2v, h2ht, h, hhv2v, htv2hhv_1, 6);
    } else if (ndiags_into == 2) {
      /* two diagonals into the far corner.
         find the face with diagonal not into the corner,
         rotate it to the right side, apply template */
      Int diag_face = -1;
      for (Int i = 0; i < 3; ++i)
        if (!diags_into[i]) diag_face = i;
      hex_bur_rot_to_right(hhv2v, diag_face);
      fill_tets_from_hex(tv2v, h2ht, h, hhv2v, htv2hhv_2, 6);
    } else {
      /* three diagonals into the far corner. */
      fill_tets_from_hex(tv2v, h2ht, h, hhv2v, htv2hhv_3, 6);
    }
  };
  parallel_for(nh, fill, "tets_from_hexes(fill)");
  return tv2v;
}

/* The symmetric algorithms place nodes at the midpoints of quads and hexes */

void tris_from_quads_symmetric(Mesh* mesh) {
  auto nold_verts = mesh->nverts();
  auto nquads = mesh->nfaces();
  auto old_coords = mesh->coords();
  auto quad_center_coords = average_field(mesh, 2, 2, old_coords);
  auto nverts = nold_verts + nquads;
  auto coords = Write<Real>(nverts * 2);
  auto ntris = nquads * 4;
  auto qv2v = mesh->ask_elem_verts();
  auto tv2v = Write<LO>(ntris * 3);
  auto f = OMEGA_H_LAMBDA(LO q) {
    auto qqv2v = gather_verts<4>(qv2v, q);
    for (Int qqt = 0; qqt < 4; ++qqt) {
      auto t = q * 4 + qqt;
      auto v0 = qqv2v[qqt];
      auto v1 = qqv2v[(qqt + 1) % 4];
      auto vq = nold_verts + q;
      tv2v[t * 3 + 0] = vq;
      tv2v[t * 3 + 1] = v0;
      tv2v[t * 3 + 2] = v1;
    }
  };
  parallel_for(nquads, f);
  map_into_range(old_coords, 0, nold_verts, coords, 2);
  map_into_range(
      quad_center_coords, nold_verts, nold_verts + nquads, coords, 2);
  Mesh new_mesh(mesh->library());
  build_from_elems_and_coords(&new_mesh, OMEGA_H_SIMPLEX, 2, tv2v, coords);
  assign(*mesh, new_mesh);
}

void tets_from_hexes_symmetric(Mesh* mesh) {
  auto nold_verts = mesh->nverts();
  auto nquads = mesh->nfaces();
  auto nhexes = mesh->nregions();
  auto old_coords = mesh->coords();
  auto quad_center_coords = average_field(mesh, 2, 3, old_coords);
  auto hex_center_coords = average_field(mesh, 3, 3, old_coords);
  auto nverts = nold_verts + nquads + nhexes;
  auto coords = Write<Real>(nverts * 3);
  auto ntets = nhexes * 6 * 4;
  auto hv2v = mesh->ask_elem_verts();
  auto hq2q = mesh->ask_down(3, 2).ab2b;
  auto tv2v = Write<LO>(ntets * 4);
  auto f = OMEGA_H_LAMBDA(LO h) {
    auto hhv2v = gather_verts<8>(hv2v, h);
    auto hhq2q = gather_down<6>(hq2q, h);
    for (Int hhq = 0; hhq < 6; ++hhq) {
      auto q = hhq2q[hhq];
      for (Int hhqt = 0; hhqt < 4; ++hhqt) {
        auto t = (h * 6 + hhq) * 4 + hhqt;
        auto v0 = hhv2v[hypercube_down_template(3, 2, hhq, hhqt)];
        auto v1 = hhv2v[hypercube_down_template(3, 2, hhq, (hhqt + 1) % 4)];
        auto vq = nold_verts + q;
        auto vh = nold_verts + nquads + h;
        tv2v[t * 4 + 0] = vq;
        tv2v[t * 4 + 1] = v1;
        tv2v[t * 4 + 2] = v0;
        tv2v[t * 4 + 3] = vh;
      }
    }
  };
  parallel_for(nhexes, f);
  map_into_range(old_coords, 0, nold_verts, coords, 3);
  map_into_range(
      quad_center_coords, nold_verts, nold_verts + nquads, coords, 3);
  map_into_range(hex_center_coords, nold_verts + nquads, nverts, coords, 3);
  Mesh new_mesh(mesh->library());
  build_from_elems_and_coords(&new_mesh, OMEGA_H_SIMPLEX, 3, tv2v, coords);
  assign(*mesh, new_mesh);
}

}  // end namespace Omega_h
