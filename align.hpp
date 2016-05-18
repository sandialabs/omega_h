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

By these definitions, an alignment code expresses
its own inverse.
*/

INLINE I8 make_code(bool is_flipped, I8 rotation, I8 which_down) {
  return static_cast<I8>((which_down << 3) | (rotation << 1) | is_flipped);
}

INLINE bool code_is_flipped(I8 code) {
  return code & 1;
}

INLINE I8 code_rotation(I8 code) {
  return (code >> 1) & 3;
}

INLINE I8 code_which_down(I8 code) {
  return (code >> 3);
}

template <Int deg>
INLINE I8 rotate_index(I8 index, I8 rotation) {
  return (index + rotation) % deg;
}

/* all the following can probably be optimized
   down to a few integer ops by an expert... */

template <Int deg>
INLINE I8 flip_index(I8 index) {
  switch(index) {
    case 1: return 2;
    case 2: return 1;
    default: return 0;
  }
}

template <Int deg>
INLINE I8 align_index(I8 index, I8 code) {
  index = rotate_index<deg>(code_rotation(code));
  if (code_is_flipped(code))
    index = flip_index<deg>(index);
  return index;
}

template <Int deg>
INLINE I8 rotation_to_first(I8 new_first) {
  return (deg - new_first) % deg;
}

/* returns the single transformation equivalent
   to applying the (code1) transformation followed
   by the (code2) one. */
template <Int deg>
INLINE I8 compound_alignments(I8 code1, I8 code2) {
  /* assuming that a code is its own inverse,
     we can look for the way to undo the compound,
     by looking at what happens to the index
     that used to be first (0) */
  I8 old_first = align_index<deg>(align_index<deg>(0, code1), code2);
  /* the inverse transformation would bring that
     index back to being the first */
  I8 rotation = rotation_to_first<deg>(old_first);
  bool is_flipped = (code_is_flipped(code1) ^ code_is_flipped(code2));
  return make_code(is_flipped, rotation, 0);
}

template <Int deg, typename T>
INLINE void rotate_adj(I8 rotation,
    T const in[], T out[]) {
  for (I8 j = 0; j < deg; ++j)
    out[rotate_index<deg>(j, rotation)] = in[j];
}

template <typename T>
INLINE void flip_adj(T adj[]) {
  swap2(adj[1], adj[2]);
}

template <Int deg, typename T>
INLINE void align_adj(I8 code,
    T const& in, T out[]) {
  rotate_adj<deg>(code_rotation(code), in, out);
  if (code_is_flipped(code))
    flip_adj<deg>(out);
}
