/* a six-bit code describes the alignment relationship
   between a simplex and a lower-dimensional simplex
   on its boundary:

   "which_down" given a canonical ordering of the
   lower-dimensional simplices on the boundary, which
   one is this one ?
   (3 bits)

The other two pieces describe alignment between two
representations of the same simplex.

  "rotation" curl-aligned, counterclockwise rotation
  (indices move up in the canonical ordering).
  0 1 2 -> 2 0 1 -> 1 2 0
  (2 bits)

  "is_flipped" only applies to triangles, swap
  the last two vertices.
  0 1 2 -> 0 2 1
  (1 bit)

We define the rotation to take place first, so a
code containing both a flip and rotation means
that to get from one entity to another one must
first rotate and then flip the vertex list.
*/

INLINE I8 make_code(bool is_flipped, Int rotation, Int which_down) {
  return static_cast<I8>((which_down << 3) | (rotation << 1) | is_flipped);
}

template <Int nverts_per_ent>
INLINE Int rotate_index(Int index, Int rotation) {
  return (index + rotation) % nverts_per_ent;
}

/* all the following can probably be optimized
   down to a few integer ops by an expert... */

INLINE Int flip_vert_index(Int index) {
  switch(index) {
    case 1: return 2;
    case 2: return 1;
    default: return 0;
  }
}

INLINE Int flip_edge_index(Int index) {
  switch(index) {
    case 0: return 2;
    case 2: return 0;
    default: return 1;
  }
}

template <Int nverts_per_ent>
INLINE Int align_vert_index(Int index, I8 code) {
  index = rotate_index<nverts_per_ent>(index, code_rotation(code));
  if (code_is_flipped(code)) {
    index = flip_vert_index(index);
  }
  return index;
}

INLINE Int align_edge_index(Int index, I8 code) {
  index = rotate_index<3>(index, code_rotation(code));
  if (code_is_flipped(code)) {
    index = flip_edge_index(index);
  }
  return index;
}

INLINE Int align_index(Int nverts_per_ent, Int index_dim, Int index, I8 code) {
  if (nverts_per_ent == 3) {
    if (index_dim == 1) {
      return align_edge_index(index, code);
    } else {
      return align_vert_index<3>(index, code);
    }
  } if (nverts_per_ent == 2) {
    return align_vert_index<2>(index, code);
  }
  NORETURN(0);
}

template <Int nverts_per_ent>
INLINE Int invert_rotation(Int rotation) {
  return (nverts_per_ent - rotation) % nverts_per_ent;
}

template <Int nverts_per_ent>
INLINE Int rotation_to_first(Int new_first) {
  return invert_rotation<nverts_per_ent>(new_first);
}

INLINE Int rotation_to_first(Int deg, Int new_first) {
  if (deg == 3) return rotation_to_first<3>(new_first);
  if (deg == 2) return rotation_to_first<2>(new_first);
  NORETURN(-1);
}

template <Int nverts_per_ent>
INLINE I8 invert_alignment(I8 code) {
  if (code_is_flipped(code))
    return code; // I think flipped codes are their own inverses
  return make_code(false,
      invert_rotation<nverts_per_ent>(code_rotation(code)), 0);
}

INLINE I8 invert_alignment(Int nverts_per_ent, I8 code) {
  if (nverts_per_ent == 3)
    return invert_alignment<3>(code);
  if (nverts_per_ent == 2)
    return invert_alignment<2>(code);
  NORETURN(0);
}

/* returns the single transformation equivalent
   to applying the (code1) transformation followed
   by the (code2) one. */
template <Int nverts_per_ent>
INLINE I8 compound_alignments(I8 code1, I8 code2) {
  /* we can look for the inverse of the compound
     by looking at what happens to the vertex
     that used to be first (0) */
  Int old_first = align_vert_index<nverts_per_ent>(
      align_vert_index<nverts_per_ent>(0, code1), code2);
  /* the inverse transformation would bring that
     vertex back to being the first */
  Int rotation = rotation_to_first<nverts_per_ent>(old_first);
  bool is_flipped = (code_is_flipped(code1) ^ code_is_flipped(code2));
  return invert_alignment<nverts_per_ent>(make_code(is_flipped, rotation, 0));
}

INLINE I8 compound_alignments(Int deg, I8 code1, I8 code2) {
  if (deg == 3) return compound_alignments<3>(code1, code2);
  if (deg == 2) return compound_alignments<2>(code1, code2);
  NORETURN(-1);
}

template <Int nverts_per_ent, typename T>
INLINE void rotate_adj(Int rotation,
    T const in[], T out[]) {
  for (I8 j = 0; j < nverts_per_ent; ++j)
    out[rotate_index<nverts_per_ent>(j, rotation)] = in[j];
}

template <typename T>
INLINE void flip_adj(T adj[]) {
  swap2(adj[1], adj[2]);
}

template <Int nverts_per_ent, typename T>
INLINE void align_adj(I8 code,
    T const in[], T out[]) {
  rotate_adj<nverts_per_ent>(code_rotation(code), in, out);
  if (code_is_flipped(code))
    flip_adj(out);
}

template <typename T>
Read<T> align_ev2v(Int deg, Read<T> ev2v, Read<I8> codes);
