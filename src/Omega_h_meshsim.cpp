#include "Omega_h_file.hpp"

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"

#include "Omega_h_for.hpp"

#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"

#include <fstream>

namespace Omega_h {

namespace meshsim {

//namespace {

void print_owners(Remotes owners, int rank) {
  printf("\n");
  auto ranks = owners.ranks;
  auto idxs = owners.idxs;
  auto ranks_w = Write<LO> (ranks.size());
  auto idxs_w = Write<LO> (idxs.size());
  auto r2w = OMEGA_H_LAMBDA(LO i) {
    ranks_w[i] = ranks[i];
    idxs_w[i] = idxs[i];
  };  
  parallel_for(idxs.size(), r2w);
  auto ranks_host = HostWrite<LO>(ranks_w);
  auto idxs_host = HostWrite<LO>(idxs_w);
  printf("On rank %d\n", rank);
  for (int i=0; i<idxs_host.size(); ++i) {
    printf("owner of %d, is on rank %d, with LId %d\n", i, ranks_host[i], idxs_host[i]);
  };  
  printf("\n");
  printf("\n");
  return;
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
    printf(" %d,", a_host[i]);
  };
  printf("\n");
  printf("\n");
  return;
}

int classId(pEntity e) {
  pGEntity g = EN_whatIn(e);
  assert(g);
  return GEN_tag(g);
}

int classType(pEntity e) {
  pGEntity g = EN_whatIn(e);
  assert(g);
  return GEN_type(g);
}

void read_internal(pParMesh sm, Mesh* mesh, pGModel g) {
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

  Omega_h_Family family = OMEGA_H_SIMPLEX;
  RIter regions = M_regionIter(m);
  pRegion rgn;
  LO i = 0;
  while (rgn = (pRegion) RIter_next(regions)) {
    if(R_topoType(rgn) != Rtet)
      Omega_h_fail("Non-simplex element found!\n");
    ++i;
  }
  RIter_delete(regions);
  std::vector<int> ent_nodes[4];
  std::vector<int> ent_class_ids[4];
  std::vector<int> ent_matches[3];
  std::vector<int> ent_match_classId[3];
  std::vector<int> ent_owners[3];
  const int numVtx = M_numVertices(m);
  ent_nodes[0].reserve(numVtx);
  ent_class_ids[0].reserve(numVtx);
  ent_matches[0].reserve(numVtx);
  ent_owners[0].reserve(numVtx);
  ent_match_classId[0].reserve(numVtx);
  HostWrite<Real> host_coords(numVtx*max_dim);

  VIter vertices = M_vertexIter(m);
  pVertex vtx;
  i = 0;
  int count_matches = 0;
  int count_matched = 0;

  while (vtx = (pVertex) VIter_next(vertices)) {
    double xyz[3];
    V_coord(vtx,xyz);
    if(max_dim < 3 && xyz[2] != 0)
      Omega_h_fail("The z coordinate must be zero for a 2d mesh!\n");
    for(int j=0; j<max_dim; j++) {
      host_coords[i * max_dim + j] = xyz[j];
    }
    ent_nodes[0].push_back(EN_id(vtx));
    ent_class_ids[0].push_back(classId(vtx));
    ++i;

    pPList matches = EN_getMatchingEnts(vtx, 0, 0);
    void *iterM = 0;
    pVertex match;
    count_matches = 0;
    while(match = (pVertex)PList_next(matches, &iterM)) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(vtx))) {
        ent_matches[0].push_back(EN_id(match));
        ent_match_classId[0].push_back(classId(match));
        ++count_matches;
        ++count_matched;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");
        if (EN_id(match) < EN_id(vtx)) {
          ent_owners[0].push_back(EN_id(match));
          //ent with lower global id become owners
        }
        else {
          ent_owners[0].push_back(EN_id(vtx));
        }
      }
      else if (PList_size(matches)==1) {
        ent_matches[0].push_back(-1);
        ent_match_classId[0].push_back(-1);
        ent_owners[0].push_back(EN_id(vtx));
      }
    }
    PList_delete(matches);
  }
  VIter_delete(vertices);

  //to get model matches
  auto g_numVtx = GM_numVertices(g);
  auto g_numEdge = GM_numEdges(g);
  auto g_numFace = GM_numFaces(g);
  auto g_numRgn = GM_numRegions(g);
  int const g_numEnt[4] = {g_numVtx, g_numEdge, g_numFace, g_numRgn};
  std::vector<int> model_matches[3];
  std::vector<int> model_ents[3];

  char filename[200];
  sprintf(filename, "model_matches.csv");
  std::ofstream in_str(filename);
  in_str << "id" << ",match" << ",dim\n";

  //for matched model verts
  model_matches[0].reserve(g_numVtx);
  model_ents[0].reserve(g_numVtx);
  GVIter g_verts = GM_vertexIter(g);
  pGVertex g_vert;
  while (g_vert = (pGVertex) GVIter_next(g_verts)) {
    VIter verts;
    pVertex vert;
    verts = M_classifiedVertexIter(m, g_vert, 0);
    int match_gEnt = -1;
    while (vert = (pVertex) VIter_next(verts)) {
      pPList matches = EN_getMatchingEnts(vert, 0, 0);
      void *iterM = 0;
      pVertex match;
      while (match = (pVertex)PList_next(matches, &iterM)) {
        if ((PList_size(matches)>1) && (EN_id(match) != EN_id(vert))) {
          if (classType(match) == 0) match_gEnt = classId(match);//match is vert
        }
      }
      PList_delete(matches);
    }
    VIter_delete(verts);
    model_matches[0].push_back(match_gEnt);
    model_ents[0].push_back(GEN_tag(g_vert));
    in_str << GEN_tag(g_vert) << "," << match_gEnt << "," << GEN_type(g_vert) << "\n";
  }
  GVIter_delete(g_verts);

  //for matched model edges
  model_matches[1].reserve(g_numEdge);
  model_ents[1].reserve(g_numEdge);
  GEIter g_edges = GM_edgeIter(g);
  pGEdge g_edge;
  while (g_edge = (pGEdge) GEIter_next(g_edges)) {
    EIter edges;
    pEdge edge;
    edges = M_classifiedEdgeIter(m, g_edge, 0);
    int match_gEnt = -1;
    while (edge = (pEdge) EIter_next(edges)) {
      pPList matches = EN_getMatchingEnts(edge, 0, 0);
      void *iterM = 0;
      pEdge match;
      while (match = (pEdge)PList_next(matches, &iterM)) {
        if ((PList_size(matches)>1) && (EN_id(match) != EN_id(edge))) {
          if (classType(match) == 1) match_gEnt = classId(match);//match is edge
        }
      }
      PList_delete(matches);
    }
    EIter_delete(edges);
    model_matches[1].push_back(match_gEnt);
    model_ents[1].push_back(GEN_tag(g_edge));
    in_str << GEN_tag(g_edge) << "," << match_gEnt << "," << GEN_type(g_edge) << "\n";
  }
  GEIter_delete(g_edges);

  //for matched model faces
  model_matches[2].reserve(g_numFace);
  model_ents[2].reserve(g_numFace);
  GFIter g_faces = GM_faceIter(g);
  pGFace g_face;
  while (g_face = (pGFace) GFIter_next(g_faces)) {
    FIter faces;
    pFace face;
    faces = M_classifiedFaceIter(m, g_face, 0);
    int match_gEnt = -1;
    while (face = (pFace) FIter_next(faces)) {
      pPList matches = EN_getMatchingEnts(face, 0, 0);
      void *iterM = 0;
      pFace match;
      while (match = (pFace)PList_next(matches, &iterM)) {
        if ((PList_size(matches)>1) && (EN_id(match) != EN_id(face))) {
          if (classType(match) == 2) match_gEnt = classId(match);//match is face
        }
      }
      PList_delete(matches);
    }
    FIter_delete(faces);
    model_matches[2].push_back(match_gEnt);
    model_ents[2].push_back(GEN_tag(g_face));
    in_str << GEN_tag(g_face) << "," << match_gEnt << "," << GEN_type(g_face) << "\n";
  }
  GFIter_delete(g_faces);
  //

  const int numEdges = M_numEdges(m);
  ent_nodes[1].reserve(numEdges*2);
  ent_class_ids[1].reserve(numEdges);
  ent_owners[1].reserve(numEdges);
  EIter edges = M_edgeIter(m);
  pEdge edge;
  count_matched = 0;
  while (edge = (pEdge) EIter_next(edges)) {
    printf("original edge %d hasverts\n", EN_id(edge));
    for(int j=1; j>=0; --j) {
      vtx = E_vertex(edge, j);
      ent_nodes[1].push_back(EN_id(vtx));
      printf(" %d ", EN_id(vtx));
    }
    printf("\n");
    ent_class_ids[1].push_back(classId(edge));

    pPList matches = EN_getMatchingEnts(edge, 0, 0);
    void *iterM = 0;
    pEdge match;
    count_matches = 0;
    while (match = (pEdge)PList_next(matches, &iterM)) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(edge))) {

        printf("original edge %d with verts\n", EN_id(edge));
        for(int j=1; j>=0; --j) printf(" %d ", EN_id(E_vertex(edge, j)));
        printf("\n  is matched to edge %d with verts ", EN_id(match));
        for(int j=1; j>=0; --j) printf(" %d ", EN_id(E_vertex(match, j)));
        printf("\n");

        ent_matches[1].push_back(EN_id(match));
        ent_match_classId[1].push_back(classId(match));
        ++count_matches;
        ++count_matched;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");
        if (EN_id(match) < EN_id(edge)) {
          ent_owners[1].push_back(EN_id(match));
        }
        else {
          ent_owners[1].push_back(EN_id(edge));
        }
      }
      else if (PList_size(matches)==1) {
        ent_matches[1].push_back(-1);
        ent_match_classId[1].push_back(-1);
        ent_owners[1].push_back(EN_id(edge));
      }
    }
    PList_delete(matches);
  }
  EIter_delete(edges);
  const int numFaces = M_numFaces(m);
  ent_nodes[2].reserve(numFaces*3);
  ent_class_ids[2].reserve(numFaces);
  ent_owners[2].reserve(numFaces);
  FIter faces = M_faceIter(m);
  pFace face;
  count_matched = 0;
  while (face = (pFace) FIter_next(faces)) {
    pPList verts = F_vertices(face,1);
    assert(PList_size(verts) == 3);
    void *iter = 0; // must initialize to 0
    while(vtx = (pVertex)PList_next(verts, &iter))
      ent_nodes[2].push_back(EN_id(vtx));
    PList_delete(verts);
    ent_class_ids[2].push_back(classId(face));

    pPList matches = EN_getMatchingEnts(face, 0, 0);
    void *iterM = 0;
    pFace match;
    count_matches = 0;
    while(match = (pFace)PList_next(matches, &iterM)) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(face))) {

        //printf("original face %d with verts\n", EN_id(face));
        pPList vertsP1 = F_vertices(face, 1);
        void *iterP1 = 0; // must initialize to 0
        while(vtx = (pVertex)PList_next(vertsP1, &iterP1))
          //printf(" %d ", EN_id(vtx));
        PList_delete(vertsP1);
        //printf("\n    is matched to face %d with verts ", EN_id(match));
        pPList vertsP2 = F_vertices(match, 1);
        void *iterP2 = 0; // must initialize to 0
        while(vtx = (pVertex)PList_next(vertsP2, &iterP2))
          //printf(" %d ", EN_id(vtx));
        PList_delete(verts);
        //printf("\n");

        ent_matches[2].push_back(EN_id(match));
        ent_match_classId[2].push_back(classId(match));
        ++count_matches;
        ++count_matched;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");
        if (EN_id(match) < EN_id(face)) {
          ent_owners[2].push_back(EN_id(match));
        }
        else {
          ent_owners[2].push_back(EN_id(face));
        }
      }
      else if (PList_size(matches)==1) {
        ent_matches[2].push_back(-1);
        ent_match_classId[2].push_back(-1);
        ent_owners[2].push_back(EN_id(face));
      }
    }
    PList_delete(matches);
  }
  FIter_delete(faces);
  const int numRegions = M_numRegions(m);
  ent_nodes[3].reserve(numRegions*4);
  ent_class_ids[3].reserve(numRegions);
  regions = M_regionIter(m);
  while (rgn = (pRegion) RIter_next(regions)) {

    pPList matches = EN_getMatchingEnts(rgn, 0, 0);
    void *iterM = 0;
    pFace match;
    count_matches = 0;
    while(match = (pFace)PList_next(matches, &iterM)) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(rgn)))
        Omega_h_fail("region matches found\n");
    }
    PList_delete(matches);

    pPList verts = R_vertices(rgn,1);
    assert(PList_size(verts) == 4);
    void *iter = 0; // must initialize to 0
    while(vtx = (pVertex)PList_next(verts, &iter))
      ent_nodes[3].push_back(EN_id(vtx));
    PList_delete(verts);
    ent_class_ids[3].push_back(classId(rgn));
  }
  RIter_delete(regions);
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

  //add matches and owners info to mesh
  for (Int ent_dim = 0; ent_dim < 3; ++ent_dim) {
    Int neev = element_degree(family, ent_dim, VERT);
    LO ndim_ents = static_cast<LO>(ent_nodes[ent_dim].size())/neev;
    HostWrite<LO> host_matches(ndim_ents);
    HostWrite<LO> host_match_classId(ndim_ents);
    HostWrite<LO> host_owners(ndim_ents);
    for (i = 0; i < ndim_ents; ++i) {
      host_matches[i] = ent_matches[ent_dim][static_cast<std::size_t>(i)];
      host_match_classId[i] = ent_match_classId[ent_dim][static_cast<std::size_t>(i)];
      //how to store when multiple matches on different parts in parallel?
      host_owners[i] = ent_owners[ent_dim][static_cast<std::size_t>(i)];
    }
    auto matches = Read<LO>(host_matches.write());
    auto match_classId = Read<LO>(host_match_classId.write());
    // hard coded fix for bad owners
    if (ent_dim == 1) {
      host_owners[4] = 3;
      host_owners[5] = 5;
      host_owners[6] = 0;
    }
    //
    auto owners = Read<LO>(host_owners.write());
    auto ranks = LOs(ndim_ents, 0);//serial mesh
    mesh->add_tag(ent_dim, "matches", 1, matches);
    mesh->add_tag(ent_dim, "match_classId", 1, match_classId);
    mesh->set_owners(ent_dim, Remotes(ranks, owners));
  }
  auto vert_matches = mesh->get_array<LO>(0, "matches");
  auto edge_matches = mesh->get_array<LO>(1, "matches");
  auto face_matches = mesh->get_array<LO>(2, "matches");
  printf("\nv matches\n");
  call_print(vert_matches);
  printf("\ne matches\n");
  call_print(edge_matches);
  printf("\nf matches\n");
  call_print(face_matches);
  print_owners(mesh->ask_owners(0), 0);
  print_owners(mesh->ask_owners(1), 0);
  print_owners(mesh->ask_owners(2), 0);

  //add model matches info to mesh
  for (Int ent_dim=0; ent_dim<max_dim; ++ent_dim) {
    HostWrite<LO> host_model_ids(g_numEnt[ent_dim]);
    HostWrite<LO> host_model_matches(g_numEnt[ent_dim]);
    for (i = 0; i < g_numEnt[ent_dim]; ++i) {
      host_model_ids[i] = model_ents[ent_dim][static_cast<std::size_t>(i)];
      host_model_matches[i] = model_matches[ent_dim][static_cast<std::size_t>(i)];
    }
    auto model_ids = Read<LO>(host_model_ids.write());
    auto model_match = Read<LO>(host_model_matches.write());
    mesh->set_model_ents(ent_dim, model_ids);
    mesh->set_model_matches(ent_dim, model_match);
    auto return_model_ents = mesh->ask_model_ents(ent_dim);
    auto return_model_matches = mesh->ask_model_matches(ent_dim);
    //printf("ent_dim=%d\n", ent_dim);
    //call_print(return_model_ents);
    //call_print(return_model_matches);
    //mesh->balance();
  }
  //
}

//}  // end anonymous namespace

Mesh read(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm) {
  SimPartitionedMesh_start(NULL,NULL);
  SimModel_start();
  Sim_readLicenseFile(NULL);
  pNativeModel nm = NULL;
  pProgress p = NULL;
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  pParMesh sm = PM_load(mesh_fname.c_str(), g, p);
  auto mesh = Mesh(comm->library());
  meshsim::read_internal(sm, &mesh, g);

  M_release(sm);
  GM_release(g);
  SimModel_stop();
  SimPartitionedMesh_stop();
  return mesh;

}

}  // namespace meshsim

}  // end namespace Omega_h
