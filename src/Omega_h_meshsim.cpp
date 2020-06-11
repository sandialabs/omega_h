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

int classId(pEntity e) {
  pGEntity g = EN_whatIn(e);
  assert(g);
  return GEN_tag(g);
}

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

  std::vector<int> elem_vertices[4];
  std::vector<int> face_vertices[2];
  std::vector<int> edge_vertices[1];
  std::vector<int> ent_class_ids[1];

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

  edge_vertices[0].reserve(numEdges*2);
  //ent_class_ids[1].reserve(numEdges);
  EIter edges = M_edgeIter(m);
  pEdge edge;
  pVertex vtx;
  int count_edge = 0;
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
  int count_tri = 0;
  int count_quad = 0;
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      count_tri += 1;
    }
    else if (F_numEdges(face) == 4) {
      count_quad += 1;
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);

  face_vertices[0].reserve(count_tri*3);
  face_vertices[1].reserve(count_quad*4);

  faces = M_faceIter(m);
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      pVertex tri_vertex;
      pPList tri_vertices = F_vertices(face,1);
      assert (PList_size(tri_vertices) == 3);
      void *iter = 0;
      while (tri_vertex = (pVertex) PList_next(tri_vertices, &iter)) {
        face_vertices[0].push_back(EN_id(tri_vertex));
      }
      PList_delete(tri_vertices);
    }
    else if (F_numEdges(face) == 4) {
      pVertex quad_vertex;
      pPList quad_vertices = F_vertices(face,1);
      assert (PList_size(quad_vertices) == 4);
      void *iter = 0;
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
  auto down = reflect_down(tri2verts, edge2vert.ab2b, vert2edge, Topo_type::triangle, Topo_type::edge);
  mesh->set_ents(Topo_type::triangle, Topo_type::edge, down);

  HostWrite<LO> host_quad2verts(count_quad*4);
  for (Int i = 0; i < count_quad; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_quad2verts[i*4 + j] =
          face_vertices[1][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto quad2verts = Read<LO>(host_quad2verts.write());
  down = reflect_down(quad2verts, edge2vert.ab2b, vert2edge, Topo_type::quadrilateral, Topo_type::edge);
  mesh->set_ents(Topo_type::quadrilateral, Topo_type::edge, down);

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
      void *iter = 0;
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[0].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rhex) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 8);
      void *iter = 0;
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[1].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rwedge) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 6);
      void *iter = 0;
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[2].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rpyramid) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 5);
      void *iter = 0;
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

  HostWrite<LO> host_tet2verts(count_tet*4);
  for (Int i = 0; i < count_tet; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_tet2verts[i*4 + j] =
          elem_vertices[0][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto tet2verts = Read<LO>(host_tet2verts.write());
  down = reflect_down(tet2verts, tri2vert.ab2b, vert2tri, Topo_type::tetrahedron, Topo_type::triangle);
  mesh->set_ents(Topo_type::tetrahedron, Topo_type::triangle, down);

  auto tet2edge = mesh->ask_down(Topo_type::tetrahedron, Topo_type::edge);
  auto tet2vtx = mesh->ask_down(Topo_type::tetrahedron, Topo_type::vertex);

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
  down = reflect_down(hex2verts, quad2vert.ab2b, vert2quad, Topo_type::hexahedron, Topo_type::quadrilateral);
  mesh->set_ents(Topo_type::hexahedron, Topo_type::quadrilateral, down);

  auto hex2edge = mesh->ask_down(Topo_type::hexahedron, Topo_type::edge);
  auto hex2vtx = mesh->ask_down(Topo_type::hexahedron, Topo_type::vertex);

  HostWrite<LO> host_wedge2verts(count_wedge*6);
  for (Int i = 0; i < count_wedge; ++i) {
    for (Int j = 0; j < 6; ++j) {
      host_wedge2verts[i*6 + j] =
          elem_vertices[2][static_cast<std::size_t>(i*6 + j)];
    }
  }
  auto wedge2verts = Read<LO>(host_wedge2verts.write());
  down = reflect_down(wedge2verts, quad2vert.ab2b, vert2quad, Topo_type::wedge, Topo_type::quadrilateral);
  mesh->set_ents(Topo_type::wedge, Topo_type::quadrilateral, down);

  auto wedge2edge = mesh->ask_down(Topo_type::wedge, Topo_type::edge);
  auto wedge2vtx = mesh->ask_down(Topo_type::wedge, Topo_type::vertex);

  down = reflect_down(wedge2verts, tri2vert.ab2b, vert2tri, Topo_type::wedge, Topo_type::triangle);
  mesh->set_ents(Topo_type::wedge, Topo_type::triangle, down);

  HostWrite<LO> host_pyramid2verts(count_pyramid*5);
  for (Int i = 0; i < count_pyramid; ++i) {
    for (Int j = 0; j < 5; ++j) {
      host_pyramid2verts[i*5 + j] =
          elem_vertices[3][static_cast<std::size_t>(i*5 + j)];
    }
  }
  auto pyramid2verts = Read<LO>(host_pyramid2verts.write());
  down = reflect_down(pyramid2verts, tri2vert.ab2b, vert2tri, Topo_type::pyramid, Topo_type::triangle);
  mesh->set_ents(Topo_type::pyramid, Topo_type::triangle, down);

  auto pyram2edge = mesh->ask_down(Topo_type::pyramid, Topo_type::edge);
  auto pyram2vtx = mesh->ask_down(Topo_type::pyramid, Topo_type::vertex);

  down = reflect_down(pyramid2verts, quad2vert.ab2b, vert2quad, Topo_type::pyramid, Topo_type::quadrilateral);
  mesh->set_ents(Topo_type::pyramid, Topo_type::quadrilateral, down);

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
  printf("wedge2vtx values from ask_down is\n");
  call_print(wedge2vtx.ab2b);
/*
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
*/
  //

/*
//ALL OK
  // test tags
  Write<LO> gravityArray(count_wedge, 10, "gravityArray");
  Read<LO> gravityArray_r(gravityArray);
  printf("ok tag1 \n");
  //mesh->add_tag<LO>(Topo_type::vertex, "gravity", 1, gravityArray_r);
  mesh->add_tag<LO>(Topo_type::wedge, "gravity", 1);
  printf("ok tag2 \n");
  mesh->set_tag<LO>(Topo_type::wedge, "gravity", gravityArray_r);
  auto test_tag = mesh->get_array<LO>(Topo_type::wedge, "gravity");
  assert(test_tag.size() == count_wedge);
  call_print(test_tag); 
  printf("ok tag3 \n");
  mesh->remove_tag(Topo_type::wedge, "gravity");
  printf("ok tag4 \n");
*/

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
  meshsim::read_internal(m, &mesh);
  M_release(m);
  GM_release(g);
  SimDiscrete_stop(0);
  SimModel_stop();
  return mesh;
}

}  // namespace meshsim

}  // end namespace Omega_h
