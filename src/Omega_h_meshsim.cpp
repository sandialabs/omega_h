#include "Omega_h_file.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_mesh.hpp"

#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_adj.hpp"
#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"

#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"

#include "SimDiscrete.h"

namespace Omega_h {

namespace meshsim {

//namespace {

/*
int classId(pEntity e) {
  pGEntity g = EN_whatIn(e);
  assert(g);
  return GEN_tag(g);
}
*/

void call_print(LOs a) {
  printf("\n");
  auto a_w = Write<LO> (a.size());
  auto r2w = OMEGA_H_LAMBDA(LO i) {
    a_w[i] = a[i];
  };
  parallel_for(a.size(), r2w);
  auto a_host = HostWrite<LO>(a_w);
  for (int i=0; i<a_host.size(); ++i) {
    printf(" %d", a_host[i]);
  };
  printf("\n");
  printf("\n");
  return;
}

void read_internal(pMesh m, Mesh* mesh) {

  (void)mesh;

  const int numVtx = M_numVertices(m);
  const int numEdges = M_numEdges(m);
  const int numFaces = M_numFaces(m);
  const int numRegions = M_numRegions(m);

/*
  Int max_dim;
  if (numRegions(m)) {
    max_dim = 3;
  } else if (numFaces(m)) {
    max_dim = 2;
  } else if (numEdges(m)) {
    max_dim = 1;
  } else {
    Omega_h_fail("There were no Elements of dimension higher than zero!\n");
  }
*/
  std::vector<int> elem_vertices[4];
  std::vector<int> face_vertices[2];
  std::vector<int> edge_vertices[1];
  //std::vector<int> down_adjs[10];
  //std::vector<int> down_codes[10];
  //std::vector<int> ent_class_ids[10];

/*
  //write vertex coords into node_coords and vertex ids into ents_nodes
  down_adjs[0].reserve(numVtx);
  ent_class_ids[0].reserve(numVtx);
  HostWrite<Real> host_coords(numVtx*max_dim);
  VIter vertices = M_vertexIter(m);
  pVertex vtx;
  i = 0;
  while ((vtx = (pVertex) VIter_next(vertices))) {
    double xyz[3];
    V_coord(vtx,xyz);
    if( max_dim < 3 && xyz[2] != 0 )
      Omega_h_fail("The z coordinate must be zero for a 2d mesh!\n");
    for(int j=0; j<max_dim; j++) {
      host_coords[i * max_dim + j] = xyz[j];
    }
    down_adjs[0].push_back(EN_id(vtx));
    ent_class_ids[0].push_back(classId(vtx));
    ++i;
  }
  VIter_delete(vertices);
*/

  //get the ids of vertices bounding each edge
  edge_vertices[0].reserve(numEdges*2);
  //ent_class_ids[1].reserve(numEdges);
  EIter edges = M_edgeIter(m);
  pEdge edge;
  pVertex vtx;
  int count_edge = 0;
  //printf(" ok1.4.2 \n");
  double xyz[3];//stores vtx coords
  while ((edge = (pEdge) EIter_next(edges))) {
    count_edge += 1;
    //printf("edge EN_id is=%d\n", EN_id(edge));
    for(int j=0; j<2; ++j) {
      vtx = E_vertex(edge,j);
      edge_vertices[0].push_back(EN_id(vtx));
      V_coord(vtx,xyz);
      //printf("vtx EN_id is=%d, x=%f, y=%f, z=%f\n", EN_id(vtx), xyz[0], xyz[1], xyz[2]);
    }
    //ent_class_ids[1].push_back(classId(edge));
  }
  EIter_delete(edges);
  //printf(" ok1.4.3 \n");
  
/*
  //Test degree outputs
  auto deg = element_degree(Topo_type::edge, Topo_type::vertex);
  printf("deg e2v=%d \n", deg);
  deg = element_degree(Topo_type::tetrahedron, Topo_type::vertex);
  printf("deg tet2v=%d \n", deg);
  deg = element_degree(Topo_type::hexahedron, Topo_type::quadrilateral);
  printf("deg hex2quad=%d \n", deg);
  deg = element_degree(Topo_type::wedge, Topo_type::triangle);
  printf("deg wed2tri=%d \n", deg);
  deg = element_degree(Topo_type::pyramid, Topo_type::edge);
  printf("deg py2ed=%d \n", deg);
  deg = element_degree(Topo_type::pyramid, Topo_type::quadrilateral);
  printf("deg py2quad=%d \n", deg);
*/
  //set Verts of mesh
  mesh->set_verts_type(numVtx);
  //
  HostWrite<LO> host_e2v(numEdges*2);
  for (Int i = 0; i < numEdges; ++i) {
    for (Int j = 0; j < 2; ++j) {
      host_e2v[i*2 + j] =
          edge_vertices[0][static_cast<std::size_t>(i*2 + j)];
    }
  }
  auto ev2v = Read<LO>(host_e2v.write());
  mesh->set_ents(Topo_type::edge, Topo_type::vertex, Adj(ev2v));

  FIter faces = M_faceIter(m);
  pFace face;
  //printf(" ok1.4.4 \n");
  int count_tri = 0;
  int count_quad = 0;
  //get ids of tris and quads from faceIDs
  //std::vector<int> face_type_ids;
  //face_type_ids.reserve(numFaces);
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      //face_type_ids[EN_id(face)] = count_tri;
      count_tri += 1;
    }
    else if (F_numEdges(face) == 4) {
      //face_type_ids[EN_id(face)] = count_quad;
      count_quad += 1;
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);

  //generate from face2verts
  face_vertices[0].reserve(count_tri*3);
  face_vertices[1].reserve(count_quad*4);

  //iterate and populate resp. vertex ids
  faces = M_faceIter(m);
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      pVertex tri_vertex;
      pPList tri_vertices = F_vertices(face,1);
      assert (PList_size(tri_vertices) == 3);
      void *iter = 0; // must initialize to 0
      while (tri_vertex = (pVertex) PList_next(tri_vertices, &iter)) {
        face_vertices[0].push_back(EN_id(tri_vertex));
      }
      PList_delete(tri_vertices);
    }
    else if (F_numEdges(face) == 4) {
      pVertex quad_vertex;
      pPList quad_vertices = F_vertices(face,1);
      assert (PList_size(quad_vertices) == 4);
      void *iter = 0; // must initialize to 0
      while (quad_vertex = (pVertex) PList_next(quad_vertices, &iter)) {
        face_vertices[1].push_back(EN_id(quad_vertex));
      }
      PList_delete(quad_vertices);
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);

  //test codes generation for tri2edge
  auto edge2vert = mesh->get_adj(Topo_type::edge, Topo_type::vertex);
  auto vert2edge = mesh->ask_up(Topo_type::vertex, Topo_type::edge);
  HostWrite<LO> host_tri2verts(count_tri*3);
  for (Int i = 0; i < count_tri; ++i) {
    for (Int j = 0; j < 3; ++j) {
      host_tri2verts[i*3 + j] =
          face_vertices[0][static_cast<std::size_t>(i*3 + j)];
    }
  }
  auto tri2verts = Read<LO>(host_tri2verts.write());
  LOs const tri_uv2v = form_uses(tri2verts, Topo_type::triangle, Topo_type::edge);
  Write<LO> te2e;
  Write<I8> tri2edg_codes;
  auto const edg2vtx_deg = element_degree(Topo_type::edge, Topo_type::vertex);
  auto const tri_2fv = get_component(tri_uv2v, edg2vtx_deg, 0);
  find_matches_ex(edg2vtx_deg, tri_2fv, tri_uv2v, edge2vert.ab2b, vert2edge, &te2e, &tri2edg_codes);
  mesh->set_ents(Topo_type::triangle, Topo_type::edge, Adj(te2e, tri2edg_codes));

  //test codes generation for quad2edge
  HostWrite<LO> host_quad2verts(count_quad*4);
  for (Int i = 0; i < count_quad; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_quad2verts[i*4 + j] =
          face_vertices[1][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto quad2verts = Read<LO>(host_quad2verts.write());
  LOs const quad_uv2v = form_uses(quad2verts, Topo_type::quadrilateral, Topo_type::edge);
  Write<LO> qe2e;
  Write<I8> quad2edg_codes;
  auto const quad_2fv = get_component(quad_uv2v, edg2vtx_deg, 0);
  find_matches_ex(edg2vtx_deg, quad_2fv, quad_uv2v, edge2vert.ab2b, vert2edge, &qe2e, &quad2edg_codes);
  mesh->set_ents(Topo_type::quadrilateral, Topo_type::edge, Adj(qe2e, quad2edg_codes));
/*
  //printf(" ok1.4.5 \n");
  //allocate memory for t2e and q2e
  down_adjs[2].reserve(count_tri*3);
  down_codes[0].reserve(count_tri*3);
  down_adjs[3].reserve(count_quad*4);
  down_codes[1].reserve(count_quad*4);
  //

  //args to generate face2edge codes
  Int which_down = 0;
  Int rotation = 0;
  bool is_flipped = false;

  //iterate and populate resp. edge ids
  faces = M_faceIter(m);
  while (face = (pFace) FIter_next(faces)) {
    //printf("face entity id=%d\n", EN_id(face));
    if (F_numEdges(face) == 3) {
      //printf ("is tri\n");
      pEdge tri_edge;
      pPList tri_edges = F_edges(face,1,0);
      assert (PList_size(tri_edges) == 3);
      void *iter = 0; // must initialize to 0
      which_down = 0;
      while (tri_edge = (pEdge) PList_next(tri_edges, &iter)) {
        down_adjs[2].push_back(EN_id(tri_edge));
        //printf("adjacent edge id=%d\n", EN_id(tri_edge));

        rotation = !F_dirUsingEdge(face, tri_edge);

        auto code = make_code(is_flipped, rotation, which_down);
        down_codes[0].push_back(code);
        ++which_down;
      }
      PList_delete(tri_edges);
    }
    else if (F_numEdges(face) == 4) {
      //printf ("is quad\n");
      pEdge quad_edge;
      pPList quad_edges = F_edges(face,1,0);
      assert (PList_size(quad_edges) == 4);
      void *iter = 0; // must initialize to 0
      which_down = 0;
      while (quad_edge = (pEdge) PList_next(quad_edges, &iter)) {
        down_adjs[3].push_back(EN_id(quad_edge));
        //printf("adjacent edge id=%d\n", EN_id(quad_edge));

        rotation = !F_dirUsingEdge(face, quad_edge);
        auto code = make_code(is_flipped, rotation, which_down);
        down_codes[1].push_back(code);
        ++which_down;
      }
      PList_delete(quad_edges);
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);
  //printf(" ok1.4.6 \n");

  HostWrite<LO> host_t2e(count_tri*3);
  HostWrite<I8> host_t2e_codes(count_tri*3);
  for (Int i = 0; i < count_tri; ++i) {
    for (Int j = 0; j < 3; ++j) {
      host_t2e[i*3 + j] =
          down_adjs[2][static_cast<std::size_t>(i*3 + j)];
      host_t2e_codes[i*3 + j] =
          down_codes[0][static_cast<std::size_t>(i*3 + j)];
    }
  }
  auto tre2e = Read<LO>(host_t2e.write());
  auto t2e_codes_name = std::string(dimensional_singular_name(Topo_type::triangle))
    + " " + dimensional_plural_name(Topo_type::edge) + " codes";
  Write<I8> t2e_codes(tre2e.size(), t2e_codes_name);
  t2e_codes = Write<I8>(host_t2e_codes);
  //mesh->set_ents(Topo_type::triangle, Topo_type::edge, Adj(tre2e, t2e_codes));

  HostWrite<LO> host_qe2e(count_quad*4);
  HostWrite<I8> host_q2e_codes(count_quad*4);
  for (Int i = 0; i < count_quad; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_qe2e[i*4 + j] =
          down_adjs[3][static_cast<std::size_t>(i*4 + j)];
      host_q2e_codes[i*4 + j] =
          down_codes[1][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto qde2e = Read<LO>(host_qe2e.write());
  auto q2e_codes_name =
std::string(dimensional_singular_name(Topo_type::quadrilateral))
    + " " + dimensional_plural_name(Topo_type::edge) + " codes";
  Write<I8> q2e_codes(qde2e.size(), q2e_codes_name);
  q2e_codes = Write<I8>(host_q2e_codes);
  //mesh->set_ents(Topo_type::quadrilateral, Topo_type::edge, Adj(qde2e, q2e_codes));
*/

  RIter regions = M_regionIter(m);
  LO count_tet = 0;
  LO count_hex = 0;
  LO count_wedge = 0;
  LO count_pyramid = 0;
  pRegion rgn;
  while (rgn = (pRegion) RIter_next(regions)) {
    if (R_topoType(rgn) == Rtet) {
      count_tet += 1;
    }
    else if (R_topoType(rgn) == Rhex) {
      count_hex += 1;
    }
    else if (R_topoType(rgn) == Rwedge) {
      count_wedge += 1;
    }
    else if (R_topoType(rgn) == Rpyramid) {
      count_pyramid += 1;
    }
    else {
      Omega_h_fail("Region is not tet, hex, wedge, or pyramid \n");
    }
  }
  RIter_delete(regions);
  printf("tet=%d, hex=%d, wedge=%d, pyramid=%d\n",
         count_tet, count_hex, count_wedge, count_pyramid);
  //

/*
  down_adjs[4].reserve(count_tet*4);
  down_adjs[5].reserve(count_hex*6);
  down_adjs[6].reserve(count_wedge*2);//tris
  down_adjs[7].reserve(count_wedge*3);
  down_adjs[8].reserve(count_pyramid*4);//tris
  down_adjs[9].reserve(count_pyramid);

  down_codes[2].reserve(count_tet*4);
  down_codes[3].reserve(count_hex*6);
  down_codes[4].reserve(count_wedge*2);//tris
  down_codes[5].reserve(count_wedge*3);
  down_codes[6].reserve(count_pyramid*4);//tris
  down_codes[7].reserve(count_pyramid);
  //printf(" ok1.4.7 \n");

  //Initialize for region2face codes
  which_down = 0;
  rotation = 0;
  is_flipped = false;
  //
  
  //iterate and populate resp. face ids
  regions = M_regionIter(m);
  while ((rgn = (pRegion) RIter_next(regions))) {
    //printf("region entity id=%d\n", EN_id(rgn));
    if (R_topoType(rgn) == Rtet) {
      //printf("is tet\n");
      pFace tri;
      pPList tris = R_faces(rgn,1);
      assert (PList_size(tris) == 4);
      void *iter = 0; // must initialize to 0
//      which_down = 0;
      while (tri = (pFace) PList_next(tris, &iter)) {
        down_adjs[4].push_back(face_type_ids[EN_id(tri)]);

        is_flipped = !R_dirUsingFace(rgn, tri);
        auto code = make_code(is_flipped, rotation, which_down);
        down_codes[2].push_back(code);
        ++which_down;

      }
      PList_delete(tris);
    }
    else if (R_topoType(rgn) == Rhex) {
      //printf("is hex\n");
      pFace quad;
      pPList quads = R_faces(rgn,1);
      assert (PList_size(quads) == 6);
      void *iter = 0; // must initialize to 0
//      which_down = 0;
//      rotation = 0;
      while (quad = (pFace) PList_next(quads, &iter)) {
        down_adjs[5].push_back(face_type_ids[EN_id(quad)]);

        is_flipped = !R_dirUsingFace(rgn, quad);
        if (which_down == 3) { rotation = 3;}//for 4 elems case, 0,1,2 not work
        else { rotation = 0;}
        auto code = make_code(is_flipped, rotation, which_down);
        //rotation = rotation_to_first(4, which_down);
        //printf ("rotation=%d\n", rotation);
        //auto code = make_code(is_flipped, rotation, 0);
        down_codes[3].push_back(code);
        ++which_down;

      }
      PList_delete(quads);
    }
    else if (R_topoType(rgn) == Rwedge) {
      //printf("is wedge\n");
      pFace w_face;
      pPList w_faces = R_faces(rgn,1);
      assert (PList_size(w_faces) == 5);
      void *iter = 0; // must initialize to 0
//      which_down = 0;
      while (w_face = (pFace) PList_next(w_faces, &iter)) {
        if (F_numEdges(w_face) == 3) { //face is tri 
          down_adjs[6].push_back(face_type_ids[EN_id(w_face)]);
	}
        else {
          //face is quad
          down_adjs[7].push_back(face_type_ids[EN_id(w_face)]);

          is_flipped = !R_dirUsingFace(rgn, w_face);
          auto code = make_code(is_flipped, rotation, which_down);
          down_codes[5].push_back(code);
          ++which_down;

	}
      }
      PList_delete(w_faces);
    }
    else if (R_topoType(rgn) == Rpyramid) {
      //printf("is pyramid\n");
      pFace p_face;
      pPList p_faces = R_faces(rgn,1);
      assert (PList_size(p_faces) == 5);
      void *iter = 0; // must initialize to 0
//      which_down = 0;
      while (p_face = (pFace) PList_next(p_faces, &iter)) {
        if (F_numEdges(p_face) == 3) { //face is tri
          down_adjs[8].push_back(face_type_ids[EN_id(p_face)]);

          is_flipped = !R_dirUsingFace(rgn, p_face);
          //rotation = rotation_to_first(3, which_down);
          //auto code = make_code(is_flipped, rotation, 0);
          auto code = make_code(is_flipped, rotation, which_down);
          down_codes[6].push_back(code);
          ++which_down;

        }
        else { // face is quad
          down_adjs[9].push_back(face_type_ids[EN_id(p_face)]);
        }
      }
      PList_delete(p_faces);
    }
    else {
      Omega_h_fail ("Region is not tet, hex, wedge, or pyramid \n");
    }
  }
  RIter_delete(regions);
  //printf(" ok1.4.8 \n");
  //

  HostWrite<LO> host_tet2tr(count_tet*4);
  //HostWrite<I8> host_tet2tr_codes(count_tet*4);
  for (Int i = 0; i < count_tet; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_tet2tr[i*4 + j] =
          down_adjs[4][static_cast<std::size_t>(i*4 + j)];
      //host_tet2tr_codes[i*4 + j] =
          //down_codes[2][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto tet2tr = Read<LO>(host_tet2tr.write());
  //auto tet2tr_codes_name =
//std::string(dimensional_singular_name(Topo_type::tetrahedron))
  //  + " " + dimensional_plural_name(Topo_type::triangle) + " codes";
  //Write<I8> tet2tr_codes(tet2tr.size(), tet2tr_codes_name);
  //tet2tr_codes = Write<I8>(host_tet2tr_codes);
  //mesh->set_ents(Topo_type::tetrahedron, Topo_type::triangle, Adj(tet2tr, tet2tr_codes));
  
  //printf(" ok1.4.8.1 \n");
  HostWrite<LO> host_hex2q(count_hex*6);
  //HostWrite<I8> host_hex2q_codes(count_hex*6);
  for (Int i = 0; i < count_hex; ++i) {
    for (Int j = 0; j < 6; ++j) {
      host_hex2q[i*6 + j] =
          down_adjs[5][static_cast<std::size_t>(i*6 + j)];
      //host_hex2q_codes[i*6 + j] =
        //  down_codes[3][static_cast<std::size_t>(i*6 + j)];
    }
  }
  auto hex2q = Read<LO>(host_hex2q.write());
  //auto hex2q_codes_name =
//std::string(dimensional_singular_name(Topo_type::hexahedron))
  //  + " " + dimensional_plural_name(Topo_type::quadrilateral) + " codes";
  //Write<I8> hex2q_codes(hex2q.size(), hex2q_codes_name);
  //hex2q_codes = Write<I8>(host_hex2q_codes);
  //mesh->set_ents(Topo_type::hexahedron, Topo_type::quadrilateral, Adj(hex2q, hex2q_codes));
  
  //printf(" ok1.4.8.2 \n");
  HostWrite<LO> host_wedge2tri(count_wedge*2);
  for (Int i = 0; i < count_wedge; ++i) {
    for (Int j = 0; j < 2; ++j) {
      host_wedge2tri[i*2 + j] =
          down_adjs[6][static_cast<std::size_t>(i*2 + j)];
    }
  }
  auto wedge2tri = Read<LO>(host_wedge2tri.write());
  //mesh->set_ents(Topo_type::wedge, Topo_type::triangle, Adj(wedge2tri));
  
  //printf(" ok1.4.8.3 \n");
  HostWrite<LO> host_wedge2quad(count_wedge*3);
  //HostWrite<I8> host_wedge2q_codes(count_wedge*3);
  for (Int i = 0; i < count_wedge; ++i) {
    for (Int j = 0; j < 3; ++j) {
      host_wedge2quad[i*3 + j] =
          down_adjs[7][static_cast<std::size_t>(i*3 + j)];
      //host_wedge2q_codes[i*3 + j] =
          //down_codes[5][static_cast<std::size_t>(i*3 + j)];
    }
  }
  auto wedge2quad = Read<LO>(host_wedge2quad.write());
  //auto wedge2q_codes_name =
//std::string(dimensional_singular_name(Topo_type::wedge))
  //  + " " + dimensional_plural_name(Topo_type::quadrilateral) + " codes";
  //Write<I8> wedge2q_codes(wedge2quad.size(), wedge2q_codes_name);
  //wedge2q_codes = Write<I8>(host_wedge2q_codes);
  //mesh->set_ents(Topo_type::wedge, Topo_type::quadrilateral, Adj(wedge2quad,
    //wedge2q_codes));
  
  //printf(" ok1.4.8.4 \n");
  HostWrite<LO> host_pyramid2tri(count_pyramid*4);
  //HostWrite<I8> host_pyram2t_codes(count_pyramid*4);
  for (Int i = 0; i < count_pyramid; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_pyramid2tri[i*4 + j] =
          down_adjs[8][static_cast<std::size_t>(i*4 + j)];
    //  host_pyram2t_codes[i*4 + j] =
      //    down_codes[6][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto pyramid2tri = Read<LO>(host_pyramid2tri.write());
  //auto pyram2t_codes_name =
//std::string(dimensional_singular_name(Topo_type::pyramid))
  //  + " " + dimensional_plural_name(Topo_type::triangle) + " codes";
  //Write<I8> pyram2t_codes(pyramid2tri.size(), pyram2t_codes_name);
  //pyram2t_codes = Write<I8>(host_pyram2t_codes);
  //mesh->set_ents(Topo_type::pyramid, Topo_type::triangle, Adj(pyramid2tri,
    //pyram2t_codes));
  
  HostWrite<LO> host_pyramid2quad(count_pyramid);
  for (Int i = 0; i < count_pyramid; ++i) {
    for (Int j = 0; j < 1; ++j) {
      host_pyramid2quad[i*1 + j] =
          down_adjs[9][static_cast<std::size_t>(i*1 + j)];
    }
  }
  auto pyramid2quad = Read<LO>(host_pyramid2quad.write());
  //mesh->set_ents(Topo_type::pyramid, Topo_type::quadrilateral, Adj(pyramid2quad));
  //printf(" ok1.4.9 \n");

  //print contents of input given to omega-> checks out
  for (std::vector<std::vector<int>>::size_type i = 1; i < 10; i++) {
    printf("from type %lu \n ", i);
    for (std::vector<int>::size_type j = 0; j < down_adjs[i].size(); j++) {
     std::cout << down_adjs[i][j] << ' ';
    }
    std::cout << std::endl;
  }
*/

  elem_vertices[0].reserve(count_tet*4);
  elem_vertices[1].reserve(count_hex*8);
  elem_vertices[2].reserve(count_wedge*6);
  elem_vertices[3].reserve(count_pyramid*5);

  regions = M_regionIter(m);
  while ((rgn = (pRegion) RIter_next(regions))) {
    if (R_topoType(rgn) == Rtet) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 4);
      void *iter = 0; // must initialize to 0
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[0].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rhex) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 8);
      void *iter = 0; // must initialize to 0
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[1].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rwedge) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 6);
      void *iter = 0; // must initialize to 0
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[2].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rpyramid) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 5);
      void *iter = 0; // must initialize to 0
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[3].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else {
      Omega_h_fail ("Region is not tet, hex, wedge, or pyramid \n");
    }
  }
  RIter_delete(regions);
  printf("elem2verts generated\n");

  printf("tri2vert values from ask_down is\n");
  auto tri2vert = mesh->ask_down(Topo_type::triangle, Topo_type::vertex);
  auto vert2tri = mesh->ask_up(Topo_type::vertex, Topo_type::triangle);

  //test codes generation for tet2tr
  HostWrite<LO> host_tet2verts(count_tet*4);
  for (Int i = 0; i < count_tet; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_tet2verts[i*4 + j] =
          elem_vertices[0][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto tet2verts = Read<LO>(host_tet2verts.write());
  LOs const tet_uv2v = form_uses(tet2verts, Topo_type::tetrahedron, Topo_type::triangle);
  Write<LO> tt2t;
  Write<I8> tet2tri_codes;
  auto const tri2vtx_deg = element_degree(Topo_type::triangle, Topo_type::vertex);
  auto const tet_2fv = get_component(tet_uv2v, tri2vtx_deg, 0);
  find_matches_ex(tri2vtx_deg, tet_2fv, tet_uv2v, tri2vert.ab2b, vert2tri, &tt2t, &tet2tri_codes);
  mesh->set_ents(Topo_type::tetrahedron, Topo_type::triangle, Adj(tt2t, tet2tri_codes));
  //mesh->set_ents(Topo_type::tetrahedron, Topo_type::triangle, Adj(tet2tr, tet2tri_codes));
  auto tet2edge = mesh->ask_down(Topo_type::tetrahedron, Topo_type::edge);
  auto tet2vtx = mesh->ask_down(Topo_type::tetrahedron, Topo_type::vertex);
  //

  //test codes generation for hex2quad
  printf("quad2vert values from ask_down is\n");
  auto quad2vert = mesh->ask_down(Topo_type::quadrilateral, Topo_type::vertex);
  auto vert2quad = mesh->ask_up(Topo_type::vertex, Topo_type::quadrilateral);

  HostWrite<LO> host_hex2verts(count_hex*8);
  for (Int i = 0; i < count_hex; ++i) {
    for (Int j = 0; j < 8; ++j) {
      host_hex2verts[i*8 + j] =
          elem_vertices[1][static_cast<std::size_t>(i*8 + j)];
    }
  }
  auto hex2verts = Read<LO>(host_hex2verts.write());
  LOs const hex_uv2v = form_uses(hex2verts, Topo_type::hexahedron, Topo_type::quadrilateral);
  Write<LO> hq2q;
  Write<I8> hex2quad_codes;
  auto const quad2vtx_deg = element_degree(Topo_type::quadrilateral, Topo_type::vertex);
  auto const hex_2fv = get_component(hex_uv2v, quad2vtx_deg, 0);
  find_matches_ex(quad2vtx_deg, hex_2fv, hex_uv2v, quad2vert.ab2b, vert2quad, &hq2q, &hex2quad_codes);
  mesh->set_ents(Topo_type::hexahedron, Topo_type::quadrilateral, Adj(hq2q, hex2quad_codes));
  //mesh->set_ents(Topo_type::hexahedron, Topo_type::quadrilateral, Adj(hex2q, hex2quad_codes));
  auto hex2edge = mesh->ask_down(Topo_type::hexahedron, Topo_type::edge);
  auto hex2vtx = mesh->ask_down(Topo_type::hexahedron, Topo_type::vertex);

  //test codes generation for wedge2quad
  HostWrite<LO> host_wedge2verts(count_wedge*6);
  for (Int i = 0; i < count_wedge; ++i) {
    for (Int j = 0; j < 6; ++j) {
      host_wedge2verts[i*6 + j] =
          elem_vertices[2][static_cast<std::size_t>(i*6 + j)];
    }
  }
  auto wedge2verts = Read<LO>(host_wedge2verts.write());
  //create wedge2quad
  LOs const wedge_uv2v = form_uses(wedge2verts, Topo_type::wedge, Topo_type::quadrilateral);
  Write<LO> wq2q;
  Write<I8> wedge2quad_codes;
  auto const wedge_2fv = get_component(wedge_uv2v, quad2vtx_deg, 0);
  find_matches_ex(quad2vtx_deg, wedge_2fv, wedge_uv2v, quad2vert.ab2b, vert2quad, &wq2q, &wedge2quad_codes);
  mesh->set_ents(Topo_type::wedge, Topo_type::quadrilateral, Adj(wq2q, wedge2quad_codes));
  //mesh->set_ents(Topo_type::wedge, Topo_type::quadrilateral, Adj(wedge2quad, wedge2quad_codes));
  auto wedge2edge = mesh->ask_down(Topo_type::wedge, Topo_type::edge);
  auto wedge2vtx = mesh->ask_down(Topo_type::wedge, Topo_type::vertex);

  //create wedge2tri
  LOs const wedgeTri_uv2v = form_uses(wedge2verts, Topo_type::wedge, Topo_type::triangle);
  Write<LO> wt2t;
  Write<I8> wedge2tri_codes;
  auto const wedgeTri_2fv = get_component(wedgeTri_uv2v, tri2vtx_deg, 0);
  find_matches_ex(tri2vtx_deg, wedgeTri_2fv, wedgeTri_uv2v, tri2vert.ab2b, vert2tri, &wt2t, &wedge2tri_codes);
  mesh->set_ents(Topo_type::wedge, Topo_type::triangle, Adj(wt2t, wedge2tri_codes));

  //test codes generation for pyramid2tri
  HostWrite<LO> host_pyramid2verts(count_pyramid*5);
  for (Int i = 0; i < count_pyramid; ++i) {
    for (Int j = 0; j < 5; ++j) {
      host_pyramid2verts[i*5 + j] =
          elem_vertices[3][static_cast<std::size_t>(i*5 + j)];
    }
  }
  auto pyramid2verts = Read<LO>(host_pyramid2verts.write());
  //create pyramid2tri
  LOs const pyramid_uv2v = form_uses(pyramid2verts, Topo_type::pyramid, Topo_type::triangle);
  Write<LO> pt2t;
  Write<I8> pyramid2tri_codes;
  auto const pyramid_2fv = get_component(pyramid_uv2v, tri2vtx_deg, 0);
  find_matches_ex(tri2vtx_deg, pyramid_2fv, pyramid_uv2v, tri2vert.ab2b, vert2tri, &pt2t, &pyramid2tri_codes);
  mesh->set_ents(Topo_type::pyramid, Topo_type::triangle, Adj(pt2t, pyramid2tri_codes));
  auto pyram2edge = mesh->ask_down(Topo_type::pyramid, Topo_type::edge);
  auto pyram2vtx = mesh->ask_down(Topo_type::pyramid, Topo_type::vertex);

  //create pyramid2quad
  LOs const pyramidQuad_uv2v = form_uses(pyramid2verts, Topo_type::pyramid, Topo_type::quadrilateral);
  Write<LO> pq2q;
  Write<I8> pyramid2quad_codes;
  auto const pyramidQuad_2fv = get_component(pyramidQuad_uv2v, quad2vtx_deg, 0);
  find_matches_ex(quad2vtx_deg, pyramidQuad_2fv, pyramidQuad_uv2v, quad2vert.ab2b, vert2quad, &pq2q, &pyramid2quad_codes);
  mesh->set_ents(Topo_type::pyramid, Topo_type::quadrilateral, Adj(pq2q, pyramid2quad_codes));

  //1 lvl dwn queries
  auto tri2edge = mesh->get_adj(Topo_type::triangle, Topo_type::edge);
  auto quad2edge = mesh->get_adj(Topo_type::quadrilateral, Topo_type::edge);
  auto tet2tri = mesh->get_adj(Topo_type::tetrahedron, Topo_type::triangle);
  auto hex2quad = mesh->get_adj(Topo_type::hexahedron, Topo_type::quadrilateral);
  auto wedge2tria = mesh->get_adj(Topo_type::wedge, Topo_type::triangle);
  auto wedge2quadr = mesh->get_adj(Topo_type::wedge, Topo_type::quadrilateral);
  auto pyramid2tria = mesh->get_adj(Topo_type::pyramid, Topo_type::triangle);
  auto pyramid2quadr = mesh->get_adj(Topo_type::pyramid, Topo_type::quadrilateral);

  //1 lvl inversion queries
  auto edge2tri = mesh->ask_up(Topo_type::edge, Topo_type::triangle);
  auto edge2quad = mesh->ask_up(Topo_type::edge, Topo_type::quadrilateral);
  auto tri2tet = mesh->ask_up(Topo_type::triangle, Topo_type::tetrahedron);
  auto quad2hex = mesh->ask_up(Topo_type::quadrilateral, Topo_type::hexahedron);
  auto tri2wedge = mesh->ask_up(Topo_type::triangle, Topo_type::wedge);
  auto quad2wedge = mesh->ask_up(Topo_type::quadrilateral, Topo_type::wedge);
  auto tri2pyramid = mesh->ask_up(Topo_type::triangle, Topo_type::pyramid);
  auto quad2pyramid = mesh->ask_up(Topo_type::quadrilateral, Topo_type::pyramid);

/*
  //inversion prints
  printf("edge2vert values returned from get_adj is\n");
  call_print(edge2vert.ab2b);
  printf("vert2edge values returned from ask_up is\n");
  call_print(vert2edge.ab2b);
  printf("vert2edge offsets returned from ask_up is\n");
  call_print(vert2edge.a2ab);
  //printf("tri2edge codes returned from get_adj is\n");
  //call_print(tri2edge.codes);
  printf("edge2tri values returned from ask_up is\n");
  call_print(edge2tri.ab2b);
  printf("edge2tri offsets returned from ask_up is\n");
  call_print(edge2tri.a2ab);
  printf("edge2quad lists returned from ask_up is\n");
  call_print(edge2quad.ab2b);
  printf("edge2quad offsets returned from ask_up is\n");
  call_print(edge2quad.a2ab);
  printf("tri2tet list returned from ask_up is\n");
  call_print(tri2tet.ab2b);
  printf("tri2tet offset returned from ask_up is\n");
  call_print(tri2tet.a2ab);
  printf("quad2hex list returned from ask_up is\n");
  call_print(quad2hex.ab2b);
  printf("quad2hex offset returned from ask_up is\n");
  call_print(quad2hex.a2ab);
  printf("tri2wedge list returned from ask_up is\n");
  call_print(tri2wedge.ab2b);
  printf("tri2wedge offset returned from ask_up is\n");
  call_print(tri2wedge.a2ab);
  printf("quad2wedge list returned from ask_up is\n");
  call_print(quad2wedge.ab2b);
  printf("quad2wedge offset returned from ask_up is\n");
  call_print(quad2wedge.a2ab);
  printf("tri2pyramid list returned from ask_up is\n");
  call_print(tri2pyramid.ab2b);
  printf("tri2pyramid offset returned from ask_up is\n");
  call_print(tri2pyramid.a2ab);
  printf("quad2pyramid list returned from ask_up is\n");
  call_print(quad2pyramid.ab2b);
  printf("quad2pyramid offset returned from ask_up is\n");
  call_print(quad2pyramid.a2ab);
*/
  //inversion assertions

  OMEGA_H_CHECK(vert2edge.ab2b.size() == 2*numEdges);
  OMEGA_H_CHECK(vert2edge.a2ab.size() == numVtx+1);
  OMEGA_H_CHECK(edge2tri.ab2b.size() == 3*count_tri);
  OMEGA_H_CHECK(edge2tri.a2ab.size() == numEdges+1);
  OMEGA_H_CHECK(edge2quad.ab2b.size() == 4*count_quad);
  OMEGA_H_CHECK(edge2quad.a2ab.size() == numEdges+1);
  OMEGA_H_CHECK(tri2tet.ab2b.size() == 4*count_tet);
  OMEGA_H_CHECK(tri2tet.a2ab.size() == count_tri+1);
  OMEGA_H_CHECK(quad2hex.ab2b.size() == count_hex*6);
  OMEGA_H_CHECK(quad2hex.a2ab.size() == count_quad+1);
  OMEGA_H_CHECK(tri2wedge.ab2b.size() == 2*count_wedge);
  OMEGA_H_CHECK(tri2wedge.a2ab.size() == count_tri+1);
  OMEGA_H_CHECK(quad2wedge.ab2b.size() == count_wedge*3);
  OMEGA_H_CHECK(quad2wedge.a2ab.size() == count_quad+1);
  OMEGA_H_CHECK(tri2pyramid.ab2b.size() == 4*count_pyramid);
  OMEGA_H_CHECK(tri2pyramid.a2ab.size() == count_tri+1);
  OMEGA_H_CHECK(quad2pyramid.ab2b.size() == count_pyramid);
  OMEGA_H_CHECK(quad2pyramid.a2ab.size() == count_quad+1);


  //transit tests
  //multi level up adj
  auto vert2tet = mesh->ask_up(Topo_type::vertex, Topo_type::tetrahedron);
  auto vert2hex = mesh->ask_up(Topo_type::vertex, Topo_type::hexahedron);
  auto vert2wedge = mesh->ask_up(Topo_type::vertex, Topo_type::wedge);
  auto vert2pyramid = mesh->ask_up(Topo_type::vertex, Topo_type::pyramid);
  auto edge2tet = mesh->ask_up(Topo_type::edge, Topo_type::tetrahedron);
  auto edge2hex = mesh->ask_up(Topo_type::edge, Topo_type::hexahedron);
  auto edge2wedge = mesh->ask_up(Topo_type::edge, Topo_type::wedge);
  auto edge2pyramid = mesh->ask_up(Topo_type::edge, Topo_type::pyramid);

  //assertions

  OMEGA_H_CHECK(tri2vert.ab2b.size() == count_tri*3);
  OMEGA_H_CHECK(quad2vert.ab2b.size() == count_quad*4);
  OMEGA_H_CHECK(tet2edge.ab2b.size() == count_tet*6);
  OMEGA_H_CHECK(tet2vtx.ab2b.size() == count_tet*4);
  OMEGA_H_CHECK(hex2edge.ab2b.size() == count_hex*12);
  OMEGA_H_CHECK(hex2vtx.ab2b.size() == count_hex*8);
  OMEGA_H_CHECK(wedge2edge.ab2b.size() == count_wedge*9);
  OMEGA_H_CHECK(wedge2vtx.ab2b.size() == count_wedge*6);
  OMEGA_H_CHECK(pyram2edge.ab2b.size() == count_pyramid*8);
  OMEGA_H_CHECK(pyram2vtx.ab2b.size() == count_pyramid*5);
  //transit prints
  printf("tri2vert values from ask_down is\n");
  call_print(tri2vert.ab2b);
  printf("quad2vert values from ask_down is\n");
  call_print(quad2vert.ab2b);
  printf("tet2edge values from ask_down is\n");
  call_print(tet2edge.ab2b);
  printf("tet2vtx values from ask_down is\n");
  call_print(tet2vtx.ab2b);
  printf("hex2edge values from ask_down is\n");
  call_print(hex2edge.ab2b);
  printf("hex2vtx values from ask_down is\n");
  call_print(hex2vtx.ab2b);
  printf("wedge2edge values from ask_down is\n");
  call_print(wedge2edge.ab2b);
  printf("wedge2vtx values from ask_down is\n");
  call_print(wedge2vtx.ab2b);
  printf("pyram2edge values from ask_down is\n");
  call_print(pyram2edge.ab2b);
  printf("pyram2vtx values from ask_down is\n");
  call_print(pyram2vtx.ab2b);
  //

/*
  //get the ids of vertices bounding each face
  ent_class_ids[2].reserve(numFaces);
  FIter faces = M_faceIter(m);
  pFace face;
  while ((face = (pFace) FIter_next(faces))) {
    ent_class_ids[2].push_back(classId(face));
  }
  FIter_delete(faces);

  //get the ids of vertices bounding each region
  ent_class_ids[3].reserve(numRegions);
  regions = M_regionIter(m);
  while ((rgn = (pRegion) RIter_next(regions))) {
    ent_class_ids[3].push_back(classId(rgn));
  }
  RIter_delete(regions);

  //flatten the down_adjs and ent_class_ids arrays
  for (Int ent_dim = max_dim; ent_dim >= 0; --ent_dim) {
    Int neev = element_degree(family, ent_dim, VERT);
    LO ndim_ents = static_cast<LO>(down_adjs[ent_dim].size()) / neev;
    HostWrite<LO> host_ev2v(ndim_ents * neev);
    HostWrite<LO> host_class_id(ndim_ents);
    for (i = 0; i < ndim_ents; ++i) {
      for (Int j = 0; j < neev; ++j) {
        host_ev2v[i * neev + j] =
            down_adjs[ent_dim][static_cast<std::size_t>(i * neev + j)];
      }
      host_class_id[i] = ent_class_ids[ent_dim][static_cast<std::size_t>(i)];
    }
    auto eqv2v = Read<LO>(host_ev2v.write());
    classify_equal_order(mesh, ent_dim, eqv2v, host_class_id.write());
  }
  finalize_classification(mesh);
*/
}

//}  // end anonymous namespace

Mesh read(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm) {
  SimPartitionedMesh_start(NULL,NULL);
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimDiscrete_start(0);
  pNativeModel nm = NULL;
  pProgress p = NULL;
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  pMesh m = M_load(mesh_fname.c_str(), g, p);
  auto mesh = Mesh(comm->library());
  meshsim::read_internal(m, &mesh);//use for fine mesh
  M_release(m);
  GM_release(g);
  SimDiscrete_stop(0);
  SimModel_stop();
  return mesh;
}

}  // namespace meshsim

}  // end namespace Omega_h
