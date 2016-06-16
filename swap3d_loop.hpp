namespace swap3d {

/* by definition, the loop vertices curl
   around the edge by the right-hand rule,
   i.e. counterclockwise when looking from
   the second edge vertex to the first. */
struct Loop {
  Int size;
  Few<LO, 2> edge_verts2verts;
  Few<LO, MAX_EDGE_SWAP> loop_verts2verts;
};

/* TODO: function to find Ring given Mesh and edge */

}//end namespace swap3d
