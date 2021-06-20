#include "Omega_h_file.hpp"

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"

#include "Omega_h_for.hpp"

#include "MeshSim.h"
#include "SimModel.h"
#include "SimUtil.h"

#include <fstream>

#include "Omega_h_array_ops.hpp"

namespace Omega_h {

namespace meshsim {

void print_matches(c_Remotes matches, int rank, int d) {
  printf("\n");
  auto root_ranks = matches.root_ranks;
  auto root_idxs = matches.root_idxs;
  auto leaf_idxs = matches.leaf_idxs;
  auto root_ranks_w = Write<LO> (root_ranks.size());
  auto root_idxs_w = Write<LO> (root_idxs.size());
  auto leaf_idxs_w = Write<LO> (leaf_idxs.size());
  auto r2w = OMEGA_H_LAMBDA(LO i) {
    root_ranks_w[i] = root_ranks[i];
    root_idxs_w[i] = root_idxs[i];
    leaf_idxs_w[i] = leaf_idxs[i];
  };  
  parallel_for(leaf_idxs.size(), r2w);
  auto root_ranks_host = HostWrite<LO>(root_ranks_w);
  auto root_idxs_host = HostWrite<LO>(root_idxs_w);
  auto leaf_idxs_host = HostWrite<LO>(leaf_idxs_w);
  printf("On rank %d, for dim = %d\n", rank, d);
  for (int i=0; i<leaf_idxs_host.size(); ++i) {
    printf("leaf id %d, root id %d, root rank %d\n", leaf_idxs_host[i],
            root_idxs_host[i], root_ranks_host[i]);
  };  
  printf("\n");
  printf("\n");
  return;
}

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

void read_matchInternal(pMesh sm, Mesh* mesh, pGModel g, CommPtr comm) {
  pMesh m = sm;
  auto rank = comm->rank();
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
  printf("read int ok-1 \n");

  Omega_h_Family family = OMEGA_H_SIMPLEX;
  RIter regions = M_regionIter(m);
  pRegion rgn;
  LO i = 0;
  while ((rgn = (pRegion) RIter_next(regions))) {
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

  //c_r
  std::vector<int> leaf_idxs[max_dim];
  std::vector<int> root_idxs[max_dim];
  std::vector<int> root_ranks[max_dim];
  //
  VIter vertices = M_vertexIter(m);
  pVertex vtx;
  i = 0;
  int count_matches = 0;
  int count_matched = 0;

  printf("read int ok0 \n");
  while ((vtx = (pVertex) VIter_next(vertices))) {
    double xyz[3];
    V_coord(vtx,xyz);
    if(max_dim < 3 && xyz[2] != 0) {
      Omega_h_fail("The z coordinate must be zero for a 2d mesh!\n");
    }
    for(int j=0; j<max_dim; j++) {
      host_coords[i * max_dim + j] = xyz[j];
    }
    ent_nodes[0].push_back(EN_id(vtx));
    ent_class_ids[0].push_back(classId(vtx));

    pPList matches = EN_getMatchingEnts(vtx, 0, 0);
    count_matches = 0;
  printf("read int ok0.0 vert %d numvtx %d rank %d\n", i, numVtx, rank);
    pVertex match;
    void *iterM = 0;
  printf("read int ok0.1 vert %d numvtx %d rank %d\n", i, numVtx, rank);
    while ((match = (pVertex) PList_next(matches, &iterM))) {
    printf("read int ok0.2 vert %d numvtx %d rank %d\n", i, numVtx, rank);
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(vtx))) {
    printf("read int ok0.3 vert %d numvtx %d rank %d\n", i, numVtx, rank);
        ent_matches[0].push_back(EN_id(match));
        ent_match_classId[0].push_back(classId(match));
    printf("read int ok0.4 vert %d numvtx %d rank %d\n", i, numVtx, rank);
        ++count_matches;
        ++count_matched;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");

        leaf_idxs[0].push_back(EN_id(vtx));
        root_idxs[0].push_back(EN_id(match));
        root_ranks[0].push_back(0);//input is a serial mesh
      }
    }
    PList_delete(matches);
    printf("read int ok0.5 vert %d numvtx %d rank %d\n", i, numVtx, rank);
    ++i;
  }
  VIter_delete(vertices);
  printf("read int ok1 \n");
  I8 mesh_IsMatched = 0;
  if (count_matched > 0) mesh_IsMatched = 1;
  mesh->set_periodic(mesh_IsMatched);

  auto g_numVtx = GM_numVertices(g);
  auto g_numEdge = GM_numEdges(g);
  auto g_numFace = GM_numFaces(g);
  auto g_numRgn = GM_numRegions(g);
  int const g_numEnt[4] = {g_numVtx, g_numEdge, g_numFace, g_numRgn};
  std::vector<int> model_matches[3];
  std::vector<int> model_ents[3];

  //for matched model verts
  model_matches[0].reserve(g_numVtx);
  model_ents[0].reserve(g_numVtx);
  GVIter g_verts = GM_vertexIter(g);
  pGVertex g_vert;
  while ((g_vert = (pGVertex) GVIter_next(g_verts))) {
    VIter verts;
    pVertex vert;
    verts = M_classifiedVertexIter(m, g_vert, 0);
    int match_gEnt = -1;
    while ((vert = (pVertex) VIter_next(verts))) {
      pPList matches = EN_getMatchingEnts(vert, 0, 0);
      void *iterM = 0;
      pVertex match;
      while ((match = (pVertex)PList_next(matches, &iterM))) {
        if ((PList_size(matches)>1) && (EN_id(match) != EN_id(vert))) {
          if (classType(match) == 0) match_gEnt = classId(match);//match is vert
        }
      }
      PList_delete(matches);
    }
    VIter_delete(verts);
    model_matches[0].push_back(match_gEnt);
    model_ents[0].push_back(GEN_tag(g_vert));
  }
  GVIter_delete(g_verts);

  //for matched model edges
  model_matches[1].reserve(g_numEdge);
  model_ents[1].reserve(g_numEdge);
  GEIter g_edges = GM_edgeIter(g);
  pGEdge g_edge;
  while ((g_edge = (pGEdge) GEIter_next(g_edges))) {
    EIter edges;
    pEdge edge;
    edges = M_classifiedEdgeIter(m, g_edge, 0);
    int match_gEnt = -1;
    while ((edge = (pEdge) EIter_next(edges))) {
      pPList matches = EN_getMatchingEnts(edge, 0, 0);
      void *iterM = 0;
      pEdge match;
      while ((match = (pEdge)PList_next(matches, &iterM))) {
        if ((PList_size(matches)>1) && (EN_id(match) != EN_id(edge))) {
          if (classType(match) == 1) match_gEnt = classId(match);//match is edge
        }
      }
      PList_delete(matches);
    }
    EIter_delete(edges);
    model_matches[1].push_back(match_gEnt);
    model_ents[1].push_back(GEN_tag(g_edge));
  }
  GEIter_delete(g_edges);

  //for matched model faces
  model_matches[2].reserve(g_numFace);
  model_ents[2].reserve(g_numFace);
  GFIter g_faces = GM_faceIter(g);
  pGFace g_face;
  while ((g_face = (pGFace) GFIter_next(g_faces))) {
    FIter faces;
    pFace face;
    faces = M_classifiedFaceIter(m, g_face, 0);
    int match_gEnt = -1;
    while ((face = (pFace) FIter_next(faces))) {
      pPList matches = EN_getMatchingEnts(face, 0, 0);
      void *iterM = 0;
      pFace match;
      while ((match = (pFace)PList_next(matches, &iterM))) {
        if ((PList_size(matches)>1) && (EN_id(match) != EN_id(face))) {
          if (classType(match) == 2) match_gEnt = classId(match);//match is face
        }
      }
      PList_delete(matches);
    }
    FIter_delete(faces);
    model_matches[2].push_back(match_gEnt);
    model_ents[2].push_back(GEN_tag(g_face));
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
  while ((edge = (pEdge) EIter_next(edges))) {
    for(int j=1; j>=0; --j) {
      vtx = E_vertex(edge, j);
      double xyz[3];
      V_coord(vtx,xyz);
      ent_nodes[1].push_back(EN_id(vtx));
    }
    ent_class_ids[1].push_back(classId(edge));

    pPList matches = EN_getMatchingEnts(edge, 0, 0);
    void *iterM = 0;
    pEdge match;
    count_matches = 0;
    while ((match = (pEdge)PList_next(matches, &iterM))) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(edge))) {
        ent_matches[1].push_back(EN_id(match));
        ent_match_classId[1].push_back(classId(match));
        ++count_matches;
        ++count_matched;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");

        leaf_idxs[1].push_back(EN_id(edge));
        root_idxs[1].push_back(EN_id(match));
        root_ranks[1].push_back(0);
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
  while ((face = (pFace) FIter_next(faces))) {
    pPList verts = F_vertices(face,1);
    assert(PList_size(verts) == 3);
    void *iter = 0; // must initialize to 0
    while ((vtx = (pVertex)PList_next(verts, &iter))) {
      ent_nodes[2].push_back(EN_id(vtx));
    }
    PList_delete(verts);
    ent_class_ids[2].push_back(classId(face));

    pPList matches = EN_getMatchingEnts(face, 0, 0);
    void *iterM = 0;
    pFace match;
    count_matches = 0;
    while ((match = (pFace)PList_next(matches, &iterM))) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(face))) {

        ent_matches[2].push_back(EN_id(match));
        ent_match_classId[2].push_back(classId(match));
        ++count_matches;
        ++count_matched;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");

        leaf_idxs[2].push_back(EN_id(face));
        root_idxs[2].push_back(EN_id(match));
        root_ranks[2].push_back(0);
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
  while ((rgn = (pRegion) RIter_next(regions))) {

    pPList matches = EN_getMatchingEnts(rgn, 0, 0);
    void *iterM = 0;
    pFace match;
    count_matches = 0;
    while ((match = (pFace)PList_next(matches, &iterM))) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(rgn))) {
        Omega_h_fail("region matches found\n");
      }
    }
    PList_delete(matches);

    pPList verts = R_vertices(rgn,1);
    assert(PList_size(verts) == 4);
    void *iter = 0; // must initialize to 0
    while ((vtx = (pVertex)PList_next(verts, &iter))) {
      ent_nodes[3].push_back(EN_id(vtx));
    }
    PList_delete(verts);
    ent_class_ids[3].push_back(classId(rgn));
  }
  RIter_delete(regions);

  auto vert_globals = Read<GO>(numVtx, 0, 1);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->set_family(family);
  mesh->set_dim(max_dim);
  mesh->set_verts(numVtx);
  mesh->add_tag(0, "global", 1, vert_globals);
  printf("read int ok2 \n");
  //add e2vrts
  for (Int ent_dim = 1; ent_dim <= max_dim; ++ent_dim) {
    Int neev = element_degree(family, ent_dim, VERT);
    LO ndim_ents = static_cast<LO>(ent_nodes[ent_dim].size()) / neev;
    HostWrite<LO> host_ev2v(ndim_ents * neev);
    for (i = 0; i < ndim_ents; ++i) {
      for (Int j = 0; j < neev; ++j) {
        host_ev2v[i * neev + j] =
          ent_nodes[ent_dim][static_cast<std::size_t>(i * neev + j)];
      }
    }
    auto ev2v = Read<LO>(host_ev2v.write());
    if (ent_dim == 1) {
      mesh->set_ents(ent_dim, Adj(ev2v));
    }
    else {
      auto ldim = ent_dim - 1;
      auto lv2v = mesh->ask_verts_of(ldim);
      auto v2l = mesh->ask_up(0, ldim);
      auto down = reflect_down(ev2v, lv2v, v2l, mesh->family(), ent_dim, ldim);
      mesh->set_ents(ent_dim, down);
    }
    printf("read int ok2 d=%d\n", ent_dim);
    mesh->add_tag(ent_dim, "global", 1, GOs(ndim_ents, 0, 1));
  }
  //
  printf("before reorder\n");
  if (!comm->reduce_and(is_sorted(vert_globals))) {
    reorder_by_globals(mesh);
  }
  printf("after reorder\n");
  mesh->add_coords(host_coords.write());
  //

  printf("read internal ok 3\n");
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
    classify_equal_order(mesh, ent_dim, host_ev2v.write(), host_class_id.write());
  }
  finalize_classification(mesh);
  printf("read internal ok 4\n");

  //add c_r input info to mesh
  for (Int ent_dim = 0; ent_dim < max_dim; ++ent_dim) {
    LO nLeaves = leaf_idxs[ent_dim].size();
    HostWrite<LO> host_leaf_idxs(nLeaves);
    HostWrite<LO> host_root_idxs(nLeaves);
    HostWrite<I32> host_root_ranks(nLeaves);
    for (i = 0; i < nLeaves; ++i) {
      host_leaf_idxs[i] = leaf_idxs[ent_dim][static_cast<std::size_t>(i)];
      host_root_idxs[i] = root_idxs[ent_dim][static_cast<std::size_t>(i)];
      host_root_ranks[i] = root_ranks[ent_dim][static_cast<std::size_t>(i)];
    }
    auto cr_leaf_idxs = Read<LO>(host_leaf_idxs.write());
    auto cr_root_idxs = Read<LO>(host_root_idxs.write());
    auto cr_root_ranks = Read<I32>(host_root_ranks.write());
    auto cr = c_Remotes(cr_leaf_idxs, cr_root_idxs, cr_root_ranks);
    mesh->set_matches(ent_dim, cr);
    auto o_cr = mesh->get_matches(ent_dim);
  }
  printf("read internal ok 5\n");

  //convert matches and owners info to device
  for (Int ent_dim = 0; ent_dim < max_dim; ++ent_dim) {
    Int neev = element_degree(family, ent_dim, VERT);
    LO ndim_ents = static_cast<LO>(ent_nodes[ent_dim].size())/neev;
    HostWrite<LO> host_matches(ndim_ents);
    HostWrite<LO> host_match_classId(ndim_ents);
    HostWrite<LO> host_owners(ndim_ents);
    for (i = 0; i < ndim_ents; ++i) {
      host_matches[i] = ent_matches[ent_dim][static_cast<std::size_t>(i)];
      host_match_classId[i] = ent_match_classId[ent_dim][static_cast<std::size_t>(i)];
      host_owners[i] = ent_owners[ent_dim][static_cast<std::size_t>(i)];
    }
    auto matches = Read<LO>(host_matches.write());
    auto match_classId = Read<LO>(host_match_classId.write());
    auto owners = Read<LO>(host_owners.write());
    auto ranks = LOs(ndim_ents, 0);//serial mesh
  }
  printf("read internal ok 6\n");

  //add model matches info to mesh and if mesh is periodic
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
  }
  printf("read internal ok 7\n");
}

void matchRead(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm, Mesh *mesh, int is_in) {
  //SimPartitionedMesh_start(NULL,NULL);
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  pNativeModel nm = NULL;
  pProgress p = NULL;
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  pMesh sm = M_load(mesh_fname.c_str(), g, p);

/*
  {
    pMesh m = PM_mesh(sm, 0);
    VIter vertices = M_vertexIter(m);
    pVertex vtx;
    int i = 0;
    while ((vtx = (pVertex) VIter_next(vertices))) {
      pPList matches = EN_getMatchingEnts(vtx, 0, 0);
      void *iterM = 0;
      pVertex match;
      while ((match = (pVertex) PList_next(matches, &iterM))) {
        printf("read int ok0.2 vert %d \n", i);
      }
      PList_delete(matches);
      ++i;
    }
    VIter_delete(vertices);
  }
*/
  if (is_in) {
    mesh->set_comm(comm);
    printf("read ok1\n");
    meshsim::read_matchInternal(sm, mesh, g, comm);
    printf("read ok2\n");
  }

  M_release(sm);
  GM_release(g);
  SimModel_stop();
  MS_exit();
  Sim_logOff();
  //SimPartitionedMesh_stop();
}

}  // namespace meshsim

}  // end namespace Omega_h
