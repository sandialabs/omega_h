#include "Omega_h_file.hpp"

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_mesh.hpp"

#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"

namespace Omega_h {

namespace meshsim {

namespace {

int classId(pEntity e) {
  pGEntity g = EN_whatIn(e);
  assert(g);
  return GEN_tag(g);
}

void read_internal(pParMesh sm, Mesh* mesh) {
  (void)mesh;
  pMesh m = PM_mesh(sm, 0);
/*
  Int max_dim;
  if (M_numRegions(m)) {
    max_dim = 3;
  } else if (M_numFaces(m)) {
    max_dim = 2;
  } else if (M_numEdges(m)) {
    max_dim = 1;
  } else {
    Omega_h_fail("There were no Elements of dimension higher than zero!\n");
  }
*/
  //get the types of elements
  //Omega_h_Family family = OMEGA_H_SIMPLEX;
  //RIter regions = M_regionIter(m);
  pRegion rgn;
  //LO i = 0;
  //while ((rgn = (pRegion) RIter_next(regions))) {
  //  if(R_topoType(rgn) != Rtet)
  //    Omega_h_fail("Non-simplex element found!\n");
  //  ++i;
  //}
  //RIter_delete(regions);
  std::vector<int> ent_nodes[10];
  //std::vector<int> ent_nodes[4];
  std::vector<int> ent_class_ids[10];
  //std::vector<int> ent_class_ids[4];
  printf(" ok1.4.1 \n");

/*
  //write vertex coords into node_coords and vertex ids into ent_nodes
  const int numVtx = M_numVertices(m);
  ent_nodes[0].reserve(numVtx);
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
    ent_nodes[0].push_back(EN_id(vtx));
    ent_class_ids[0].push_back(classId(vtx));
    ++i;
  }
  VIter_delete(vertices);
*/
  //get the ids of vertices bounding each edge
  const int numEdges = M_numEdges(m);
  ent_nodes[1].reserve(numEdges*2);
  ent_class_ids[1].reserve(numEdges);
  EIter edges = M_edgeIter(m);
  pEdge edge;
  pVertex vtx;
  printf(" ok1.4.2 \n");
  while ((edge = (pEdge) EIter_next(edges))) {
    for(int j=1; j>=0; --j) {
      vtx = E_vertex(edge,j);
      ent_nodes[1].push_back(EN_id(vtx));
    }
    ent_class_ids[1].push_back(classId(edge));
  }
  EIter_delete(edges);
  printf(" ok1.4.3 \n");

  //pass vectors to set_ents
  HostWrite<LO> host_e2v(numEdges*2);
  for (Int i = 0; i < numEdges; ++i) {
    for (Int j = 0; j < 2; ++j) {
      host_e2v[i*2 + j] =
          ent_nodes[1][static_cast<std::size_t>(i*2 + j)];
    }
  }
  auto e2v = Read<LO>(host_e2v.write()); //This is LOs
  //Mesh::set_ents(1, 0, e2v);
  //may use the Adj(down) constructor to pass adj graph

  //get the ids of edges bounding each triangle
  //get the ids of edges bounding each quadrilateral
  FIter faces = M_faceIter(m);
  pFace face;
  printf(" ok1.4.4 \n");
  const int numFaces = M_numFaces(m);
  //alloc for ids of tris
  std::vector<int> tri_ids[1];
  tri_ids[0].reserve(numFaces);
  //alloc for ids of quads
  std::vector<int> quad_ids[1];
  quad_ids[0].reserve(numFaces);
  //count no. of tris and quads
  int count_tri = 0;
  int count_quad = 0;
  //get ids of tris and quads from faceIDs & total
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      tri_ids[0][EN_id(face)] = count_tri;
      //increment for next tri
      count_tri += 1;
    }
    else if (F_numEdges(face) == 4) {
      quad_ids[0][EN_id(face)] = count_quad;
      //increment for next quad
      count_quad += 1;
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);
  printf(" ok1.4.5 \n");
  //
  //allocate memory for t2e and q2e
  ent_nodes[2].reserve(count_tri*3);
  ent_nodes[3].reserve(count_quad*4);
  //
  
  //iterate and populate resp. edge ids
  faces = M_faceIter(m);
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      pEdge tri_edge;
      pPList tri_edges = F_edges(face,1,0);
      assert (PList_size(tri_edges) == 3);
      void *iter = 0; // must initialize to 0
      while (tri_edge = (pEdge) PList_next(tri_edges, &iter))
        ent_nodes[2].push_back(EN_id(tri_edge));
      PList_delete(tri_edges);
    }
    else if (F_numEdges(face) == 4) {
      pEdge quad_edge;
      pPList quad_edges = F_edges(face,1,0);
      assert (PList_size(quad_edges) == 4);
      void *iter = 0; // must initialize to 0
      //check if PList_next or E...
      while (quad_edge = (pEdge) PList_next(quad_edges, &iter))
        ent_nodes[3].push_back(EN_id(quad_edge));
      PList_delete(quad_edges);
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);
  //set ents at end
  printf(" ok1.4.6 \n");

  //pass vectors to set_ents
  HostWrite<LO> host_t2e(count_tri*3);
  for (Int i = 0; i < count_tri; ++i) {
    for (Int j = 0; j < 3; ++j) {
      host_t2e[i*3 + j] =
          ent_nodes[2][static_cast<std::size_t>(i*3 + j)];
    }
  }
  auto t2e = Read<LO>(host_t2e.write()); //This is LOs
  //Mesh::set_ents(2, 1, t2e);

  HostWrite<LO> host_q2e(count_quad*4);
  for (Int i = 0; i < count_quad; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_q2e[i*4 + j] =
          ent_nodes[3][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto q2e = Read<LO>(host_q2e.write()); //This is LOs
  //Mesh::set_ents(3, 1, q2e);

  //get the ids of tris bounding each tet
  //get the ids of quads bounding each hex
  //get the ids of tris bounding each wedge
  //get the ids of quads bounding each wedge
  //get the ids of tris bounding each pyramid
  //get the ids of quads bounding each pyramid
  RIter regions = M_regionIter(m);
  //count tets, hexs, wedges and pyramids
  LO count_tet = 0;
  LO count_hex = 0;
  LO count_wedge = 0;
  LO count_pyramid = 0;
  while (rgn = (pRegion) RIter_next(regions)) {
    if (R_topoType(rgn) == Rtet) {
      count_tet += 1;
    }
    // check which exact keyword for other types
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
  //

  printf(" ok1.4.6 \n");
  //allocate memory for t2t. h2q, w2t, w2q, p2t, p2q
  ent_nodes[4].reserve(count_tet*4);
  ent_nodes[5].reserve(count_hex*6);
  ent_nodes[6].reserve(count_wedge*2);//tris
  ent_nodes[7].reserve(count_wedge*3);
  ent_nodes[8].reserve(count_pyramid*4);//tris
  ent_nodes[9].reserve(count_pyramid);
  //
  printf(" ok1.4.7 \n");
  
  //iterate and populate resp. face ids
  regions = M_regionIter(m);
  while ((rgn = (pRegion) RIter_next(regions))) {
    //Tets
    if (R_topoType(rgn) == Rtet) {
      printf(" tet \n");
      pFace tri;
      pPList tris = R_faces(rgn,1);
      assert (PList_size(tris) == 4);
      void *iter = 0; // must initialize to 0
      while (tri = (pFace) PList_next(tris, &iter))
        ent_nodes[4].push_back(tri_ids[0][EN_id(tri)]);
        //ent_nodes[4].push_back(EN_id(tri));
      PList_delete(tris);
    }
    //Hexs
    else if (R_topoType(rgn) == Rhex) {
      printf(" hex \n");
      pFace quad;
      pPList quads = R_faces(rgn,1); 
      assert (PList_size(quads) == 6);
      void *iter = 0; // must initialize to 0
      while (quad = (pFace) PList_next(quads, &iter))
        ent_nodes[5].push_back(quad_ids[0][EN_id(quad)]);
        //ent_nodes[5].push_back(EN_id(quad));
      PList_delete(quads);
    }
    //Wedges
    else if (R_topoType(rgn) == Rwedge) {
      printf(" wedge \n");
      pFace w_face;
      pPList w_faces = R_faces(rgn,1);
      assert (PList_size(w_faces) == 5);
      void *iter = 0; // must initialize to 0
      while (w_face = (pFace) PList_next(w_faces, &iter)) {
        if (F_numEdges(w_face) == 3) { //face is tri 
          ent_nodes[6].push_back(tri_ids[0][EN_id(w_face)]);
	  //ent_nodes[6].push_back(EN_id(w_face));
	}
        else { //face is quad
          ent_nodes[7].push_back(quad_ids[0][EN_id(w_face)]);
	  //ent_nodes[7].push_back(EN_id(w_face));
	}
      }
      PList_delete(w_faces);
    }
    //Pyramids
    else if (R_topoType(rgn) == Rpyramid) {
      printf(" pyramid 1\n");
      pFace p_face;
      pPList p_faces = R_faces(rgn,1);
      assert (PList_size(p_faces) == 5);
      void *iter = 0; // must initialize to 0
      while (p_face = (pFace) PList_next(p_faces, &iter)) {
        if (F_numEdges(p_face) == 3) { //face is tri
          printf(" pyramid,tri 4.t \n");
          ent_nodes[8].push_back(tri_ids[0][EN_id(p_face)]);
          //ent_nodes[8].push_back(EN_id(p_face));
        }
        else { // face is quad
          printf(" pyramid,quad 4.q, count_py=%d \n", count_pyramid);
          ent_nodes[9].push_back(quad_ids[0][EN_id(p_face)]);
          //ent_nodes[9].push_back(EN_id(p_face));
        }
      }
      PList_delete(p_faces);
    }
    else {
      Omega_h_fail ("Region is not tet, hex, wedge, or pyramid \n");
    }
  printf(" ok1.4.7.1 \n");
  }
  RIter_delete(regions);
  printf(" ok1.4.8 \n");

  //
  //pass vectors to set_ents
  HostWrite<LO> host_tet2tr(count_tet*4);
  for (Int i = 0; i < count_tet; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_tet2tr[i*4 + j] =
          ent_nodes[4][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto tet2tr = Read<LO>(host_tet2tr.write()); //This is LOs
  //Mesh::set_ents(4, 2, tet2tr);
  
  HostWrite<LO> host_hex2q(count_hex*6);
  for (Int i = 0; i < count_hex; ++i) {
    for (Int j = 0; j < 6; ++j) {
      host_tet2tr[i*6 + j] =
          ent_nodes[5][static_cast<std::size_t>(i*6 + j)];
    }
  }
  auto hex2q = Read<LO>(host_hex2q.write()); //This is LOs
  //Mesh::set_ents(5, 3, hex2q);
  
  HostWrite<LO> host_wedge2tri(count_wedge*2);
  for (Int i = 0; i < count_wedge; ++i) {
    for (Int j = 0; j < 2; ++j) {
      host_wedge2tri[i*2 + j] =
          ent_nodes[6][static_cast<std::size_t>(i*2 + j)];
    }
  }
  auto wedge2tri = Read<LO>(host_wedge2tri.write()); //This is LOs
  //Mesh::set_ents(6, 2, wedge2tri);
  
  HostWrite<LO> host_wedge2quad(count_wedge*3);
  for (Int i = 0; i < count_wedge; ++i) {
    for (Int j = 0; j < 3; ++j) {
      host_wedge2quad[i*3 + j] =
          ent_nodes[7][static_cast<std::size_t>(i*3 + j)];
    }
  }
  auto wedge2quad = Read<LO>(host_wedge2quad.write()); //This is LOs
  //Mesh::set_ents(6, 3, wedge2quad);
  
  HostWrite<LO> host_pyramid2tri(count_pyramid*4);
  for (Int i = 0; i < count_pyramid; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_pyramid2tri[i*4 + j] =
          ent_nodes[8][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto pyramid2tri = Read<LO>(host_pyramid2tri.write()); //This is LOs
  //Mesh::set_ents(7, 2, pyramid2tri);
  
  HostWrite<LO> host_pyramid2quad(count_pyramid);
  for (Int i = 0; i < count_pyramid; ++i) {
    for (Int j = 0; j < 1; ++j) {
      host_pyramid2quad[i*1 + j] =
          ent_nodes[9][static_cast<std::size_t>(i*1 + j)];
    }
  }
  auto pyramid2quad = Read<LO>(host_pyramid2quad.write()); //This is LOs
  //Mesh::set_ents(8, 3, pyramid2quad);
/*

  const int numFaces = M_numFaces(m);
  ent_nodes[2].reserve(numFaces*3);
  ent_class_ids[2].reserve(numFaces);

  while ((face = (pFace) FIter_next(faces))) {
    pPList verts = F_vertices(face,1);
    assert(PList_size(verts) == 3);
    void *iter = 0; // must initialize to 0
    while((vtx = (pVertex)PList_next(verts, &iter)))
      ent_nodes[2].push_back(EN_id(vtx));
    PList_delete(verts);
    ent_class_ids[2].push_back(classId(face));
  }
  FIter_delete(faces);

  //get the ids of vertices bounding each face
  const int numFaces = M_numFaces(m);
  ent_nodes[2].reserve(numFaces*3);
  ent_class_ids[2].reserve(numFaces);
  FIter faces = M_faceIter(m);
  pFace face;
  while ((face = (pFace) FIter_next(faces))) {
    pPList verts = F_vertices(face,1);
    assert(PList_size(verts) == 3);
    void *iter = 0; // must initialize to 0
    while((vtx = (pVertex)PList_next(verts, &iter)))
      ent_nodes[2].push_back(EN_id(vtx));
    PList_delete(verts);
    ent_class_ids[2].push_back(classId(face));
  }
  FIter_delete(faces);

  //get the ids of vertices bounding each region
  const int numRegions = M_numRegions(m);
  ent_nodes[3].reserve(numRegions*4);
  ent_class_ids[3].reserve(numRegions);
  regions = M_regionIter(m);
  while ((rgn = (pRegion) RIter_next(regions))) {
    pPList verts = R_vertices(rgn,1);
    assert(PList_size(verts) == 4);
    void *iter = 0; // must initialize to 0
    while((vtx = (pVertex)PList_next(verts, &iter)))
      ent_nodes[3].push_back(EN_id(vtx));
    PList_delete(verts);
    ent_class_ids[3].push_back(classId(rgn));
  }
  RIter_delete(regions);

  //flatten the ent_nodes and ent_class_ids arrays
  for (Int ent_dim = max_dim; ent_dim >= 0; --ent_dim) {
    Int neev = element_degree(family, ent_dim, VERT);
    LO ndim_ents = static_cast<LO>(ent_nodes[ent_dim].size()) / neev;
    HostWrite<LO> host_ev2v(ndim_ents * neev);
    HostWrite<LO> host_class_id(ndim_ents);
    for (i = 0; i < ndim_ents; ++i) {
      for (Int j = 0; j < neev; ++j) {
        host_ev2v[i * neev + j] =
            ent_nodes[ent_dim][static_cast<std::size_t>(i * neev + j)];
      }
      host_class_id[i] = ent_class_ids[ent_dim][static_cast<std::size_t>(i)];
    }
    auto eqv2v = Read<LO>(host_ev2v.write());
    if (ent_dim == max_dim) {
      build_from_elems_and_coords(
          mesh, family, max_dim, eqv2v, host_coords.write());
    }
    classify_equal_order(mesh, ent_dim, eqv2v, host_class_id.write());
  }
  finalize_classification(mesh);
*/
}

}  // end anonymous namespace

Mesh read(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm) {
  SimPartitionedMesh_start(NULL,NULL);
  SimModel_start();
  Sim_readLicenseFile(NULL);
  pNativeModel nm = NULL;
  pProgress p = NULL;
  printf(" ok1.1 \n");
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  printf(" ok1.2 \n");
  pParMesh sm = PM_load(mesh_fname.c_str(), g, p);
  printf(" ok1.3 \n");
  auto mesh = Mesh(comm->library());
  printf(" ok1.4 \n");
  meshsim::read_internal(sm, &mesh);
  printf(" ok1.5 \n");
  M_release(sm);
  GM_release(g);
  SimModel_stop();
  SimPartitionedMesh_stop();
  return mesh;
}

}  // namespace meshsim

}  // end namespace Omega_h
