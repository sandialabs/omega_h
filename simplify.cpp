namespace simplify {

/* The algorithms here are described in:

   Dompierre, Julien, et al.
   "How to Subdivide Pyramids, Prisms, and Hexahedra into Tetrahedra."
   Proceedings of the 8th International Meshing Roundtable. 1999.
*/

namespace {

INLINE Int find_min(LO const v[], Int n) {
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

INLINE bool lt(Int v_i, Int v_j, Int v_k, Int v_l) {
  return min2(v_i, v_j) < min2(v_k, v_l);
}

INLINE void rot(LO v[], Int n) {
  LO tmp = v[n - 1];
  for (Int i = 0; i < n - 1; ++i) v[i + 1] = v[i];
  v[0] = tmp;
}

INLINE void rot_ntimes(LO v[], Int nv, Int ntimes) {
  for (Int i = 0; i < ntimes; ++i) rot(v, nv);
}

INLINE void rot_to_first(LO v[], Int nv, Int first) {
  rot_ntimes(v, nv, ((nv - first) % nv));
}

/* quad to tri template */
CONSTANT Int const qtv2qqv[2][3] = {{0, 1, 2}, {2, 3, 0}};

/* below are the hex-to-tet templates for the
   four unique cases identified by Dompierre et al.
   They list multiple templates for each case,
   here we'll be a bit lazy and just pick one each */

/* tets from a hex with no diagonals
   into the back-upper-right corner */
CONSTANT Int const htv2hhv_0[5][4] = {
    {0, 1, 2, 5}, {0, 2, 7, 5}, {0, 2, 3, 7}, {0, 5, 7, 4}, {2, 7, 5, 6}};
/* tets from a hex with 1 diagonal
   into the back-upper-right corner,
   on the right face */
CONSTANT Int const htv2hhv_1[6][4] = {{0, 5, 7, 4}, {0, 1, 7, 5}, {1, 6, 7, 5},
                                      {0, 7, 2, 3}, {0, 7, 1, 2}, {1, 7, 6, 2}};
/* tets from a hex with 2 diagonals
   into the back-upper-right corner,
   none on the right face */
CONSTANT Int const htv2hhv_2[6][4] = {{0, 4, 5, 6}, {0, 3, 7, 6}, {0, 7, 4, 6},
                                      {0, 1, 2, 5}, {0, 3, 6, 2}, {0, 6, 5, 2}};
/* tets from a hex with 3 diagonals
   into the back-upper-right corner */
CONSTANT Int const htv2hhv_3[6][4] = {{0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6},
                                      {0, 5, 6, 4}, {1, 5, 6, 0}, {1, 6, 2, 0}};

CONSTANT Int const hex_flip_pairs[4][2] = {{0, 4}, {3, 5}, {1, 7}, {2, 6}};

DEVICE void flip_hex(LO hhv2v[]) {
  for (Int i = 0; i < 4; ++i)
    swap2(hhv2v[hex_flip_pairs[i][0]], hhv2v[hex_flip_pairs[i][1]]);
}

/* faces adjacent to the back-upper-right
   corner, starting with the right face
   and curling around the centroidal XYZ axis.
   also, numbered with the corner vertex first */
CONSTANT Int const hex_bur_faces[3][4] = {
    {6, 5, 1, 2}, {6, 2, 3, 7}, {6, 7, 4, 5}};

/* the vertices that rotate amonst one another
   when a hex is rotated around its centroidal
   XYZ axis (the line between the front-lower-left
   corner and the back-upper-right corner). */
CONSTANT Int const hex_bur_ring[6] = {1, 2, 3, 7, 5, 4};

DEVICE void hex_bur_rot_ntimes(LO hhv2v[], Int ntimes) {
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
DEVICE void hex_bur_rot_to_right(LO hhv2v[], Int new_right) {
  hex_bur_rot_ntimes(hhv2v, ((3 - new_right) % 3));
}
}

LOs tris_from_quads(LOs qv2v) {
  CHECK(qv2v.size() % 4 == 0);
  LO nq = qv2v.size() / 4;
  LO nt = nq * 2;
  Write<LO> tv2v(nt * 3);
  auto f = LAMBDA(LO q) {
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
  parallel_for(nq, f);
  return tv2v;
}

DEVICE void tets_from_hex_1(LO h, LOs hv2v, LO hhv2v[], Int diags_into[],
                            Int& ndiags_into) {
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
     back-upper-left corner */
  for (Int i = 0; i < 3; ++i) {
    diags_into[i] = lt(hhv2v[hex_bur_faces[i][0]], hhv2v[hex_bur_faces[i][2]],
                       hhv2v[hex_bur_faces[i][1]], hhv2v[hex_bur_faces[i][3]]);
  }
  ndiags_into = diags_into[0] + diags_into[1] + diags_into[2];
}

DEVICE void fill_tets_from_hex(Write<LO> tv2v, LOs h2ht, LO h, LO const hhv2v[],
                               Int const case_template[][4], Int nhht) {
  LO t = h2ht[h];
  for (Int i = 0; i < nhht; ++i) {
    for (Int j = 0; j < 4; ++j) {
      tv2v[t * 4 + j] = hhv2v[case_template[i][j]];
    }
    ++t;
  }
}

LOs tets_from_hexes(LOs hv2v) {
  CHECK(hv2v.size() % 8 == 0);
  LO nh = hv2v.size() / 8;
  Write<LO> degrees(nh);
  auto count = LAMBDA(LO h) {
    LO hhv2v[8];
    Int diags_into[3];
    Int ndiags_into;
    tets_from_hex_1(h, hv2v, hhv2v, diags_into, ndiags_into);
    if (ndiags_into == 0)
      degrees[h] = 5;
    else
      degrees[h] = 6;
  };
  parallel_for(nh, count);
  auto h2ht = offset_scan(LOs(degrees));
  LO nt = h2ht.last();
  Write<LO> tv2v(nt * 4);
  auto fill = LAMBDA(LO h) {
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
  parallel_for(nh, fill);
  return tv2v;
}

}  // end namespace simplify
