#include "Omega_h_file.hpp"

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"

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
  //get the types of elements
  //Omega_h_Family family = OMEGA_H_SIMPLEX;
  //RIter regions = M_regionIter(m);
  pRegion rgn;
  LO i = 0;
  //while ((rgn = (pRegion) RIter_next(regions))) {
  //  if(R_topoType(rgn) != Rtet)
  //    Omega_h_fail("Non-simplex element found!\n");
  //  ++i;
  //}
  //RIter_delete(regions);
  std::vector<int> ent_nodes[9];
  //std::vector<int> ent_nodes[4];
  std::vector<int> ent_class_ids[9];
  //std::vector<int> ent_class_ids[4];

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
  while ((edge = (pEdge) EIter_next(edges))) {
    for(int j=1; j>=0; --j) {
      vtx = E_vertex(edge,j);
      ent_nodes[1].push_back(EN_id(vtx));
    }
    ent_class_ids[1].push_back(classId(edge));
  }
  EIter_delete(edges);

  //below code may have syntactic/keyword issues
  //get the ids of edges bounding each triangle
  //get the ids of edges bounding each quadrilateral
  FIter faces = M_faceIter(m);
  pFace face;
  //count tris and quads
  LO count_tri = 0;
  LO count_quad = 0;
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
  //

  //allocate memory for t2e and q2e
  ent_nodes[2].reserve(count_tri*3);
  ent_nodes[3].reserve(count_quad*4);
  //
  
  //iterate and populate resp. edge ids
  FIter faces = M_faceIter(m);
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      pPList tri_edges = F_edges(face,1);
      assert (PList_size(tri_edges) == 3);
      void *iter = 0; // must initialize to 0
      while (tri_edge = (pEdge) EList_next(tri_edges, &iter))
        ent_nodes[2].push_back(EN_id(tri_edge));
      PList_delete(tri_edges);
    }
    else if (F_numEdges(face) == 4) {
      pPList quad_edges = F_edges(face,1);
      assert (PList_size(quad_edges) == 4);
      void *iter = 0; // must initialize to 0
      //check if PList_next or E...
      while (quad_edge = (pEdge) EList_next(quad_edges, &iter))
        ent_nodes[3].push_back(EN_id(quad_edge));
      PList_delete(quad_edges);
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);
  //set ents at end

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

  //allocate memory for t2t. h2q, w2t, w2q, p2t, p2q
  ent_nodes[4].reserve(count_tet*4);
  ent_nodes[5].reserve(count_hex*6);
  ent_nodes[6].reserve(count_wedge*2);//tris
  ent_nodes[7].reserve(count_wedge*3);
  ent_nodes[8].reserve(count_pyramid*4);//tris
  ent_nodes[9].reserve(count_pyramid);
  //
  
  //iterate and populate resp. face ids
  RIter regions = M_regionIter(m);
  while ((rgn = (pRegion) RIter_next(regions))) {
    //Tets
    if (R_topoType(rgn) == Rtet) {
      pPList tris = R_faces(rgn,1);
      assert (PList_size(tris) == 4);
      void *iter = 0; // must initialize to 0
      while (tri = (pFace) FList_next(tris, &iter))
        ent_nodes[4].push_back(EN_id(tri));
      PList_delete(tris);
    }
    //Hexs
    else if (R_topoType(rgn) == Rhex) {
      pPList quads = R_faces(rgn,1); 
      assert (PList_size(quads) == 6);
      void *iter = 0; // must initialize to 0
      while (quad = (pFace) FList_next(quads, &iter))
        ent_nodes[5].push_back(EN_id(quad));
      PList_delete(quads);
    }
    //Wedges
    else if (R_topoType(rgn) == Rwedge) {
      pPList faces = R_faces(rgn,1);
      assert (PList_size(faces) == 5);
      void *iter = 0; // must initialize to 0
      while (face = (pFace) FList_next(faces, &iter)) {
        if (F_numEdges(face) == 3)
          { ent_nodes[6].push_back(EN_id(face)) };
        else
          { ent_nodes[7].push_back(EN_id(face)) };
      }
      PList_delete(faces);
    }
    //Pyramids
    else if (R_topoType(rgn) == Rpyramid) {
      pPList faces = R_faces(rgn,1);
      assert (PList_size(faces) == 5);
      void *iter = 0; // must initialize to 0
      while (face = (pFace) FList_next(faces, &iter)) {
        if (F_numEdges(face) == 3)
          ent_nodes[8].push_back(EN_id(face));
        else
          ent_nodes[9].push_back(EN_id(face));
      }
      PList_delete(faces);
    }
    else {
      Omega_h_fail ("Region is not tet, hex, wedge, or pyramid \n");
    }
  }
  RIter_delete(regions);
  //
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
  printf("ok1");
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  printf("ok2");
  pParMesh sm = PM_load(mesh_fname.c_str(), g, p);
  printf("ok3");
  auto mesh = Mesh(comm->library());
  printf("ok4");
  meshsim::read_internal(sm, &mesh);
  printf("ok5");
  M_release(sm);
  GM_release(g);
  SimModel_stop();
  SimPartitionedMesh_stop();
  return mesh;
}

}  // namespace meshsim

}  // end namespace Omega_h
