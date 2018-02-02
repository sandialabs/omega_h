#ifndef SWAP3D_TABLES_HPP
#define SWAP3D_TABLES_HPP

namespace Omega_h {

namespace swap3d {

enum {
  MAX_EDGE_SWAP = 7,
  MAX_UNIQUE_TRIS = 35,
  MAX_UNIQUE_EDGES = 26,
  MAX_CONFIGS = 42
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

OMEGA_H_CONSTANT_DATA static Int const swap_mesh_sizes[MAX_EDGE_SWAP + 1] = {
    0  // 0
    ,
    0  // 1
    ,
    0  // 2
    ,
    1  // 3
    ,
    2  // 4
    ,
    3  // 5
    ,
    4  // 6
    ,
    5  // 7
};

OMEGA_H_CONSTANT_DATA static Int const swap_mesh_counts[MAX_EDGE_SWAP + 1] = {
    0  // 0
    ,
    0  // 1
    ,
    0  // 2
    ,
    1  // 3
    ,
    2  // 4
    ,
    5  // 5
    ,
    14  // 6
    ,
    42  // 7
};

OMEGA_H_CONSTANT_DATA static Int const triangles_3[1][3] = {{0, 1, 2}};

OMEGA_H_CONSTANT_DATA static Int const meshes_3[1] = {0};

OMEGA_H_CONSTANT_DATA static Int const triangles_4[4][3] = {
    {0, 1, 2}, {0, 2, 3}, {0, 1, 3}, {1, 2, 3}};

OMEGA_H_CONSTANT_DATA static Int const meshes_4[2 * 2] = {0, 1, 2, 3};

OMEGA_H_CONSTANT_DATA static Int const triangles_5[10][3] = {{0, 1, 2},
    {0, 2, 3}, {0, 3, 4}, {0, 1, 4}, {1, 3, 4}, {1, 2, 3}, {2, 3, 4}, {0, 2, 4},
    {0, 1, 3}, {1, 2, 4}};

OMEGA_H_CONSTANT_DATA static Int const meshes_5[5 * 3] = {
    0, 1, 2, 3, 4, 5, 0, 6, 7, 2, 5, 8, 3, 6, 9};

OMEGA_H_CONSTANT_DATA static Int const triangles_6[20][3] = {{0, 1, 2},
    {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 2, 5}, {2, 4, 5}, {2, 3, 4}, {0, 3, 5},
    {3, 4, 5}, {0, 2, 4}, {2, 3, 5}, {1, 2, 3}, {0, 1, 3}, {0, 1, 5}, {1, 4, 5},
    {1, 3, 4}, {0, 1, 4}, {1, 3, 5}, {1, 2, 4}, {1, 2, 5}};

OMEGA_H_CONSTANT_DATA static Int const meshes_6[14 * 4] = {0, 1, 2, 3, 0, 4, 5,
    6, 0, 1, 7, 8, 0, 3, 6, 9, 0, 4, 8, 10, 2, 3, 11, 12, 11, 13, 14, 15, 7, 8,
    11, 12, 3, 11, 15, 16, 8, 11, 13, 17, 6, 13, 14, 18, 3, 6, 16, 18, 5, 6, 13,
    19, 8, 10, 13, 19};

OMEGA_H_CONSTANT_DATA static Int const triangles_7[35][3] = {{0, 1, 2},
    {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 5, 6}, {0, 3, 6}, {3, 5, 6}, {3, 4, 5},
    {0, 4, 6}, {4, 5, 6}, {0, 3, 5}, {3, 4, 6}, {0, 2, 4}, {2, 3, 4}, {0, 2, 6},
    {2, 5, 6}, {2, 4, 5}, {0, 2, 5}, {2, 4, 6}, {2, 3, 5}, {2, 3, 6}, {0, 1, 3},
    {1, 2, 3}, {0, 1, 4}, {1, 3, 4}, {0, 1, 6}, {1, 5, 6}, {1, 4, 5}, {0, 1, 5},
    {1, 4, 6}, {1, 3, 5}, {1, 3, 6}, {1, 2, 4}, {1, 2, 5}, {1, 2, 6}};

OMEGA_H_CONSTANT_DATA static Int const meshes_7[42 * 5] = {0, 1, 2, 3, 4, 0, 1,
    5, 6, 7, 0, 1, 2, 8, 9, 0, 1, 4, 7, 10, 0, 1, 5, 9, 11, 0, 3, 4, 12, 13, 0,
    13, 14, 15, 16, 0, 8, 9, 12, 13, 0, 4, 13, 16, 17, 0, 9, 13, 14, 18, 0, 7,
    14, 15, 19, 0, 4, 7, 17, 19, 0, 6, 7, 14, 20, 0, 9, 11, 14, 20, 2, 3, 4, 21,
    22, 5, 6, 7, 21, 22, 2, 8, 9, 21, 22, 4, 7, 10, 21, 22, 5, 9, 11, 21, 22, 3,
    4, 22, 23, 24, 22, 24, 25, 26, 27, 8, 9, 22, 23, 24, 4, 22, 24, 27, 28, 9,
    22, 24, 25, 29, 7, 22, 25, 26, 30, 4, 7, 22, 28, 30, 6, 7, 22, 25, 31, 9,
    11, 22, 25, 31, 3, 4, 13, 23, 32, 13, 25, 26, 27, 32, 8, 9, 13, 23, 32, 4,
    13, 27, 28, 32, 9, 13, 25, 29, 32, 13, 16, 25, 26, 33, 4, 13, 16, 28, 33,
    13, 15, 16, 25, 34, 9, 13, 18, 25, 34, 7, 19, 25, 26, 33, 4, 7, 19, 28, 33,
    7, 15, 19, 25, 34, 6, 7, 20, 25, 34, 9, 11, 20, 25, 34};

/* array [8] of pointer to array [3] of Int const */
OMEGA_H_CONSTANT_DATA static swap_tri_t const* const
    swap_triangles[MAX_EDGE_SWAP + 1] = {nullptr, nullptr, nullptr, triangles_3,
        triangles_4, triangles_5, triangles_6, triangles_7};

OMEGA_H_CONSTANT_DATA static Int const* const swap_meshes[MAX_EDGE_SWAP + 1] = {
    nullptr, nullptr, nullptr, meshes_3, meshes_4, meshes_5, meshes_6,
    meshes_7};

typedef int IntPair[2];

OMEGA_H_CONSTANT_DATA static IntPair const unique_edges_4[2] = {{2, 0}, {1, 3}};

OMEGA_H_CONSTANT_DATA static IntPair const unique_edges_5[7] = {
    {2, 0}, {3, 0}, {1, 4}, {1, 3}, {4, 2}, {0, 3}, {3, 1}};

OMEGA_H_CONSTANT_DATA static IntPair const unique_edges_6[15] = {{2, 0}, {3, 0},
    {4, 0}, {2, 5}, {2, 4}, {3, 5}, {0, 4}, {4, 2}, {5, 3}, {0, 3}, {3, 1},
    {1, 5}, {1, 4}, {4, 1}, {5, 2}};

OMEGA_H_CONSTANT_DATA static IntPair const unique_edges_7[26] = {{2, 0}, {3, 0},
    {4, 0}, {5, 0}, {3, 6}, {3, 5}, {4, 6}, {0, 5}, {5, 3}, {6, 4}, {0, 4},
    {2, 4}, {4, 2}, {2, 6}, {2, 5}, {5, 2}, {6, 3}, {0, 3}, {1, 3}, {3, 1},
    {1, 4}, {4, 1}, {1, 6}, {1, 5}, {5, 1}, {6, 2}};

OMEGA_H_CONSTANT_DATA static IntPair const* const
    unique_edges[MAX_EDGE_SWAP + 1] = {nullptr, nullptr, nullptr, nullptr,
        unique_edges_4, unique_edges_5, unique_edges_6, unique_edges_7};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_4_0[1] = {0};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_4_1[1] = {1};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_5_0[2] = {0, 1};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_5_1[2] = {2, 3};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_5_2[2] = {0, 4};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_5_3[2] = {5, 6};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_5_4[2] = {2, 4};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_0[3] = {0, 1, 2};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_1[3] = {0, 3, 4};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_2[3] = {0, 1, 5};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_3[3] = {0, 6, 7};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_4[3] = {0, 3, 8};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_5[3] = {9, 2, 10};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_6[3] = {10, 11, 12};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_7[3] = {9, 5, 10};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_8[3] = {6, 10, 13};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_9[3] = {8, 10, 11};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_10[3] = {7, 11, 12};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_11[3] = {6, 7, 12};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_12[3] = {4, 14, 11};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_6_13[3] = {8, 14, 11};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_0[4] = {0, 1, 2, 3};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_1[4] = {0, 1, 4, 5};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_2[4] = {0, 1, 2, 6};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_3[4] = {0, 1, 7, 8};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_4[4] = {0, 1, 4, 9};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_5[4] = {0, 10, 3, 11};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_6[4] = {0, 12, 13, 14};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_7[4] = {0, 10, 6, 11};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_8[4] = {0, 7, 12, 15};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_9[4] = {0, 9, 12, 13};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_10[4] = {0, 8, 13, 14};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_11[4] = {0, 7, 8, 14};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_12[4] = {0, 5, 16, 13};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_13[4] = {0, 9, 16, 13};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_14[4] = {17, 2, 3, 18};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_15[4] = {17, 4, 5, 18};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_16[4] = {17, 2, 6, 18};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_17[4] = {7, 8, 17, 18};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_18[4] = {17, 4, 9, 18};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_19[4] = {10, 3, 19, 20};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_20[4] = {19, 21, 22, 23};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_21[4] = {10, 6, 19, 20};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_22[4] = {7, 19, 21, 24};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_23[4] = {9, 19, 21, 22};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_24[4] = {8, 19, 22, 23};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_25[4] = {7, 8, 19, 23};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_26[4] = {5, 16, 19, 22};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_27[4] = {9, 16, 19, 22};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_28[4] = {10, 3, 12, 20};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_29[4] = {12, 22, 23, 20};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_30[4] = {10, 6, 12, 20};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_31[4] = {7, 12, 20, 24};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_32[4] = {9, 12, 22, 20};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_33[4] = {12, 15, 22, 23};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_34[4] = {7, 12, 15, 23};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_35[4] = {12, 14, 25, 22};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_36[4] = {9, 12, 25, 22};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_37[4] = {8, 15, 22, 23};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_38[4] = {7, 8, 15, 23};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_39[4] = {8, 14, 25, 22};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_40[4] = {5, 16, 25, 22};

OMEGA_H_CONSTANT_DATA static Int const edges2unique_7_41[4] = {9, 16, 25, 22};

OMEGA_H_CONSTANT_DATA static Int const* const edges2unique_4[2] = {
    edges2unique_4_0, edges2unique_4_1};

OMEGA_H_CONSTANT_DATA static Int const* const edges2unique_5[5] = {
    edges2unique_5_0, edges2unique_5_1, edges2unique_5_2, edges2unique_5_3,
    edges2unique_5_4};

OMEGA_H_CONSTANT_DATA static Int const* const edges2unique_6[14] = {
    edges2unique_6_0, edges2unique_6_1, edges2unique_6_2, edges2unique_6_3,
    edges2unique_6_4, edges2unique_6_5, edges2unique_6_6, edges2unique_6_7,
    edges2unique_6_8, edges2unique_6_9, edges2unique_6_10, edges2unique_6_11,
    edges2unique_6_12, edges2unique_6_13};

OMEGA_H_CONSTANT_DATA static Int const* const edges2unique_7[42] = {
    edges2unique_7_0, edges2unique_7_1, edges2unique_7_2, edges2unique_7_3,
    edges2unique_7_4, edges2unique_7_5, edges2unique_7_6, edges2unique_7_7,
    edges2unique_7_8, edges2unique_7_9, edges2unique_7_10, edges2unique_7_11,
    edges2unique_7_12, edges2unique_7_13, edges2unique_7_14, edges2unique_7_15,
    edges2unique_7_16, edges2unique_7_17, edges2unique_7_18, edges2unique_7_19,
    edges2unique_7_20, edges2unique_7_21, edges2unique_7_22, edges2unique_7_23,
    edges2unique_7_24, edges2unique_7_25, edges2unique_7_26, edges2unique_7_27,
    edges2unique_7_28, edges2unique_7_29, edges2unique_7_30, edges2unique_7_31,
    edges2unique_7_32, edges2unique_7_33, edges2unique_7_34, edges2unique_7_35,
    edges2unique_7_36, edges2unique_7_37, edges2unique_7_38, edges2unique_7_39,
    edges2unique_7_40, edges2unique_7_41};

OMEGA_H_CONSTANT_DATA static Int const* const* const
    edges2unique[MAX_EDGE_SWAP + 1] = {nullptr, nullptr, nullptr, nullptr,
        edges2unique_4, edges2unique_5, edges2unique_6, edges2unique_7};

OMEGA_H_CONSTANT_DATA static Int const nedges[MAX_EDGE_SWAP + 1] = {
    0, 0, 0, 0, 1, 2, 3, 4};

}  // end namespace swap3d

}  // end namespace Omega_h

#endif
