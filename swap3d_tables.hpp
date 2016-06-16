namespace swap3d {

enum {
  MAX_EDGE_SWAP = 7,
  MAX_UNIQUE_TRIS = 35
};

typedef Int const swap_tri_t[3];

/* tables hardcoding all possible 2D meshes (triangulations) of
   an N-sided polygon, where N ranges from 3 to 7.

   for a given polygon, many possible triangulations may share
   the same triangle, so the list of unique triangles is stored
   separately.
   this is useful because it allows caching of quality metrics
   computed based on the unique triangle which can then be
   reused when evaluating the quality of triangulations. */

CONSTANT static Int const swap_mesh_sizes[MAX_EDGE_SWAP+1] =
{0 //0
,0 //1
,0 //2
,1 //3
,2 //4
,3 //5
,4 //6
,5 //7
};

CONSTANT static Int const swap_mesh_counts[MAX_EDGE_SWAP+1] =
{0  //0
,0  //1
,0  //2
,1  //3
,2  //4
,5  //5
,14 //6
,42 //7
};

CONSTANT static Int const triangles_3[1][3] = {{0,1,2}};

CONSTANT static Int const meshes_3[1] =
{0
};

CONSTANT static Int const triangles_4[4][3] =
{{0,1,2}
,{0,2,3}
,{0,1,3}
,{1,2,3}
};

CONSTANT static Int const meshes_4[2 * 2] =
{0,1
,2,3
};

CONSTANT static Int const triangles_5[10][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,1,4}
,{1,3,4}
,{1,2,3}
,{2,3,4}
,{0,2,4}
,{0,1,3}
,{1,2,4}
};

CONSTANT static Int const meshes_5[5 * 3] =
{0,1,2
,3,4,5
,0,6,7
,2,5,8
,3,6,9
};

CONSTANT static Int const triangles_6[20][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,4,5}
,{0,2,5}
,{2,4,5}
,{2,3,4}
,{0,3,5}
,{3,4,5}
,{0,2,4}
,{2,3,5}
,{1,2,3}
,{0,1,3}
,{0,1,5}
,{1,4,5}
,{1,3,4}
,{0,1,4}
,{1,3,5}
,{1,2,4}
,{1,2,5}
};

CONSTANT static Int const meshes_6[14 * 4] =
{0, 1 ,2 ,3
,0, 4 ,5 ,6
,0, 1 ,7 ,8
,0, 3 ,6 ,9
,0, 4 ,8 ,10
,2, 3 ,11,12
,11,13,14,15
,7 ,8 ,11,12
,3 ,11,15,16
,8 ,11,13,17
,6 ,13,14,18
,3 ,6 ,16,18
,5 ,6 ,13,19
,8 ,10,13,19
};

CONSTANT static Int const triangles_7[35][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,4,5}
,{0,5,6}
,{0,3,6}
,{3,5,6}
,{3,4,5}
,{0,4,6}
,{4,5,6}
,{0,3,5}
,{3,4,6}
,{0,2,4}
,{2,3,4}
,{0,2,6}
,{2,5,6}
,{2,4,5}
,{0,2,5}
,{2,4,6}
,{2,3,5}
,{2,3,6}
,{0,1,3}
,{1,2,3}
,{0,1,4}
,{1,3,4}
,{0,1,6}
,{1,5,6}
,{1,4,5}
,{0,1,5}
,{1,4,6}
,{1,3,5}
,{1,3,6}
,{1,2,4}
,{1,2,5}
,{1,2,6}
};

CONSTANT static Int const meshes_7[42 * 5] =
{0 ,1 ,2 ,3 ,4
,0 ,1 ,5 ,6 ,7
,0 ,1 ,2 ,8 ,9
,0 ,1 ,4 ,7 ,10
,0 ,1 ,5 ,9 ,11
,0 ,3 ,4 ,12,13
,0 ,13,14,15,16
,0 ,8 ,9 ,12,13
,0 ,4 ,13,16,17
,0 ,9 ,13,14,18
,0 ,7 ,14,15,19
,0 ,4 ,7 ,17,19
,0 ,6 ,7 ,14,20
,0 ,9 ,11,14,20
,2 ,3 ,4 ,21,22
,5 ,6 ,7 ,21,22
,2 ,8 ,9 ,21,22
,4 ,7 ,10,21,22
,5 ,9 ,11,21,22
,3 ,4 ,22,23,24
,22,24,25,26,27
,8 ,9 ,22,23,24
,4 ,22,24,27,28
,9 ,22,24,25,29
,7 ,22,25,26,30
,4 ,7 ,22,28,30
,6 ,7 ,22,25,31
,9 ,11,22,25,31
,3 ,4 ,13,23,32
,13,25,26,27,32
,8 ,9 ,13,23,32
,4 ,13,27,28,32
,9 ,13,25,29,32
,13,16,25,26,33
,4 ,13,16,28,33
,13,15,16,25,34
,9 ,13,18,25,34
,7 ,19,25,26,33
,4 ,7 ,19,28,33
,7 ,15,19,25,34
,6 ,7 ,20,25,34
,9 ,11,20,25,34
};

/* array [8] of pointer to array [3] of Int const */
CONSTANT static swap_tri_t const* const swap_triangles[MAX_EDGE_SWAP+1] =
{0
,0
,0
,triangles_3
,triangles_4
,triangles_5
,triangles_6
,triangles_7
};

CONSTANT static Int const* const swap_meshes[MAX_EDGE_SWAP+1] =
{0
,0
,0
,meshes_3
,meshes_4
,meshes_5
,meshes_6
,meshes_7
};

/* the following tables were auto-generated long ago
   based on the tables above.
   don't modify manually !

   they describe, for each possible triangulation
   of an N-sided polygon, the interior edges of that
   triangulation.

   these will be used for full topology edge swapping
   to create intermediate entities in the interior of the
   cavity.

   we don't bother identifying unique edges because no caching
   is necessary for these intermediate entities */

CONSTANT static Int const edges_4_0[1 * 2] =
{2,0
};

CONSTANT static Int const edges_4_1[1 * 2] =
{1,3
};

CONSTANT static Int const edges_5_0[2 * 2] =
{2,0
,3,0
};

CONSTANT static Int const edges_5_1[2 * 2] =
{1,4
,1,3
};

CONSTANT static Int const edges_5_2[2 * 2] =
{2,0
,4,2
};

CONSTANT static Int const edges_5_3[2 * 2] =
{0,3
,3,1
};

CONSTANT static Int const edges_5_4[2 * 2] =
{1,4
,4,2
};

CONSTANT static Int const edges_6_0[3 * 2] =
{2,0
,3,0
,4,0
};

CONSTANT static Int const edges_6_1[3 * 2] =
{2,0
,2,5
,2,4
};

CONSTANT static Int const edges_6_2[3 * 2] =
{2,0
,3,0
,3,5
};

CONSTANT static Int const edges_6_3[3 * 2] =
{2,0
,0,4
,4,2
};

CONSTANT static Int const edges_6_4[3 * 2] =
{2,0
,2,5
,5,3
};

CONSTANT static Int const edges_6_5[3 * 2] =
{0,3
,4,0
,3,1
};

CONSTANT static Int const edges_6_6[3 * 2] =
{3,1
,1,5
,1,4
};

CONSTANT static Int const edges_6_7[3 * 2] =
{0,3
,3,5
,3,1
};

CONSTANT static Int const edges_6_8[3 * 2] =
{0,4
,3,1
,4,1
};

CONSTANT static Int const edges_6_9[3 * 2] =
{5,3
,3,1
,1,5
};

CONSTANT static Int const edges_6_10[3 * 2] =
{4,2
,1,5
,1,4
};

CONSTANT static Int const edges_6_11[3 * 2] =
{0,4
,4,2
,1,4
};

CONSTANT static Int const edges_6_12[3 * 2] =
{2,4
,5,2
,1,5
};

CONSTANT static Int const edges_6_13[3 * 2] =
{5,3
,5,2
,1,5
};

CONSTANT static Int const edges_7_0[4 * 2] =
{2,0
,3,0
,4,0
,5,0
};

CONSTANT static Int const edges_7_1[4 * 2] =
{2,0
,3,0
,3,6
,3,5
};

CONSTANT static Int const edges_7_2[4 * 2] =
{2,0
,3,0
,4,0
,4,6
};

CONSTANT static Int const edges_7_3[4 * 2] =
{2,0
,3,0
,0,5
,5,3
};

CONSTANT static Int const edges_7_4[4 * 2] =
{2,0
,3,0
,3,6
,6,4
};

CONSTANT static Int const edges_7_5[4 * 2] =
{2,0
,0,4
,5,0
,2,4
};

CONSTANT static Int const edges_7_6[4 * 2] =
{2,0
,4,2
,2,6
,2,5
};

CONSTANT static Int const edges_7_7[4 * 2] =
{2,0
,0,4
,4,6
,2,4
};

CONSTANT static Int const edges_7_8[4 * 2] =
{2,0
,0,5
,4,2
,5,2
};

CONSTANT static Int const edges_7_9[4 * 2] =
{2,0
,6,4
,4,2
,2,6
};

CONSTANT static Int const edges_7_10[4 * 2] =
{2,0
,5,3
,2,6
,2,5
};

CONSTANT static Int const edges_7_11[4 * 2] =
{2,0
,0,5
,5,3
,2,5
};

CONSTANT static Int const edges_7_12[4 * 2] =
{2,0
,3,5
,6,3
,2,6
};

CONSTANT static Int const edges_7_13[4 * 2] =
{2,0
,6,4
,6,3
,2,6
};

CONSTANT static Int const edges_7_14[4 * 2] =
{0,3
,4,0
,5,0
,1,3
};

CONSTANT static Int const edges_7_15[4 * 2] =
{0,3
,3,6
,3,5
,1,3
};

CONSTANT static Int const edges_7_16[4 * 2] =
{0,3
,4,0
,4,6
,1,3
};

CONSTANT static Int const edges_7_17[4 * 2] =
{0,5
,5,3
,0,3
,1,3
};

CONSTANT static Int const edges_7_18[4 * 2] =
{0,3
,3,6
,6,4
,1,3
};

CONSTANT static Int const edges_7_19[4 * 2] =
{0,4
,5,0
,3,1
,1,4
};

CONSTANT static Int const edges_7_20[4 * 2] =
{3,1
,4,1
,1,6
,1,5
};

CONSTANT static Int const edges_7_21[4 * 2] =
{0,4
,4,6
,3,1
,1,4
};

CONSTANT static Int const edges_7_22[4 * 2] =
{0,5
,3,1
,4,1
,5,1
};

CONSTANT static Int const edges_7_23[4 * 2] =
{6,4
,3,1
,4,1
,1,6
};

CONSTANT static Int const edges_7_24[4 * 2] =
{5,3
,3,1
,1,6
,1,5
};

CONSTANT static Int const edges_7_25[4 * 2] =
{0,5
,5,3
,3,1
,1,5
};

CONSTANT static Int const edges_7_26[4 * 2] =
{3,5
,6,3
,3,1
,1,6
};

CONSTANT static Int const edges_7_27[4 * 2] =
{6,4
,6,3
,3,1
,1,6
};

CONSTANT static Int const edges_7_28[4 * 2] =
{0,4
,5,0
,4,2
,1,4
};

CONSTANT static Int const edges_7_29[4 * 2] =
{4,2
,1,6
,1,5
,1,4
};

CONSTANT static Int const edges_7_30[4 * 2] =
{0,4
,4,6
,4,2
,1,4
};

CONSTANT static Int const edges_7_31[4 * 2] =
{0,5
,4,2
,1,4
,5,1
};

CONSTANT static Int const edges_7_32[4 * 2] =
{6,4
,4,2
,1,6
,1,4
};

CONSTANT static Int const edges_7_33[4 * 2] =
{4,2
,5,2
,1,6
,1,5
};

CONSTANT static Int const edges_7_34[4 * 2] =
{0,5
,4,2
,5,2
,1,5
};

CONSTANT static Int const edges_7_35[4 * 2] =
{4,2
,2,5
,6,2
,1,6
};

CONSTANT static Int const edges_7_36[4 * 2] =
{6,4
,4,2
,6,2
,1,6
};

CONSTANT static Int const edges_7_37[4 * 2] =
{5,3
,5,2
,1,6
,1,5
};

CONSTANT static Int const edges_7_38[4 * 2] =
{0,5
,5,3
,5,2
,1,5
};

CONSTANT static Int const edges_7_39[4 * 2] =
{5,3
,2,5
,6,2
,1,6
};

CONSTANT static Int const edges_7_40[4 * 2] =
{3,5
,6,3
,6,2
,1,6
};

CONSTANT static Int const edges_7_41[4 * 2] =
{6,4
,6,3
,6,2
,1,6
};

CONSTANT static Int const* const edges_4[2] =
{edges_4_0
,edges_4_1
};

CONSTANT static Int const* const edges_5[5] =
{edges_5_0
,edges_5_1
,edges_5_2
,edges_5_3
,edges_5_4
};

CONSTANT static Int const* const edges_6[14] =
{edges_6_0
,edges_6_1
,edges_6_2
,edges_6_3
,edges_6_4
,edges_6_5
,edges_6_6
,edges_6_7
,edges_6_8
,edges_6_9
,edges_6_10
,edges_6_11
,edges_6_12
,edges_6_13
};

CONSTANT static Int const* const edges_7[42] =
{edges_7_0
,edges_7_1
,edges_7_2
,edges_7_3
,edges_7_4
,edges_7_5
,edges_7_6
,edges_7_7
,edges_7_8
,edges_7_9
,edges_7_10
,edges_7_11
,edges_7_12
,edges_7_13
,edges_7_14
,edges_7_15
,edges_7_16
,edges_7_17
,edges_7_18
,edges_7_19
,edges_7_20
,edges_7_21
,edges_7_22
,edges_7_23
,edges_7_24
,edges_7_25
,edges_7_26
,edges_7_27
,edges_7_28
,edges_7_29
,edges_7_30
,edges_7_31
,edges_7_32
,edges_7_33
,edges_7_34
,edges_7_35
,edges_7_36
,edges_7_37
,edges_7_38
,edges_7_39
,edges_7_40
,edges_7_41
};

CONSTANT static Int const* const* const swap_int_edges[MAX_EDGE_SWAP + 1] =
{0
,0
,0
,0
,edges_4
,edges_5
,edges_6
,edges_7
};

CONSTANT static Int const swap_nint_edges[MAX_EDGE_SWAP + 1] =
{0
,0
,0
,0
,1
,2
,3
,4
};

} //end namespace swap3d
