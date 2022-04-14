#include "Omega_h_file.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_array_ops.hpp"

#include "Omega_h_simConfig.h"
#include "MeshSim.h"
#include "SimModel.h"
#include "SimUtil.h"

#include <fstream>

namespace {
  int classId(pEntity e) {
    pGEntity g = EN_whatIn(e);
    assert(g);
    return GEN_tag(g);
  }

  int classType(pEntity e) {
    pGEntity g = EN_whatIn(e);
    assert(g);
    assert((0 <= GEN_type(g)) && (3 >= GEN_type(g)));
    return GEN_type(g);
  }
}

namespace Omega_h {

namespace meshsim {

void read_matchInternal(pMesh sm, Mesh* mesh, pGModel g, CommPtr comm) {
  pMesh m = sm;
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
  while ((rgn = (pRegion) RIter_next(regions))) {
    if(R_topoType(rgn) != Rtet)
      Omega_h_fail("Non-simplex element found!\n");
    ++i;
  }
  RIter_delete(regions);
  std::vector<int> ent_nodes[4];
  std::vector<int> ent_class_ids[4];
  const int numVtx = M_numVertices(m);
  ent_nodes[0].reserve(numVtx);
  ent_class_ids[0].reserve(numVtx);
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

#if (SIMMODSUITE_MAJOR_VERSION <= 14 && SIMMODSUITE_MINOR_VERSION < 190921) || \
    (SIMMODSUITE_MAJOR_VERSION >= 16 && SIMMODSUITE_MINOR_VERSION >= 210623)
    pPList matches = EN_getMatchingEnts(vtx, NULL);
#else
    pPList matches = EN_getMatchingEnts(vtx, NULL, 0);
#endif
    count_matches = 0;
    pVertex match;
    void *iterM = 0;
    if (matches) { // this check fixes the segfault as
                   // EN_getMatchingEnts returns null, not an empty list, when
                   // there are no matches
      assert(PList_size(matches) > 1);
      while ((match = (pVertex) PList_next(matches, &iterM))) {
        if (EN_id(match) != EN_id(vtx)) {
          ++count_matches;
          ++count_matched;
          if (count_matches > 1) {
            Omega_h_fail("Error:matches per entity > 1\n");
          }

          leaf_idxs[0].push_back(EN_id(vtx));
          root_idxs[0].push_back(EN_id(match));
          root_ranks[0].push_back(0);//input is a serial mesh
        }
      }
    }
    PList_delete(matches);
    ++i;
  }
  VIter_delete(vertices);
  if (count_matched == 0) {
    Omega_h_fail("Error: Input mesh is not matched\n");
  }

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
#if (SIMMODSUITE_MAJOR_VERSION <= 14 && SIMMODSUITE_MINOR_VERSION < 190921) || \
    (SIMMODSUITE_MAJOR_VERSION >= 16 && SIMMODSUITE_MINOR_VERSION >= 210623)
      pPList matches = EN_getMatchingEnts(vert, NULL);
#else
      pPList matches = EN_getMatchingEnts(vert, NULL, 0);
#endif
      void *iterM = 0;
      pVertex match;
      if (matches) {
        assert(PList_size(matches) > 1);
        while ((match = (pVertex)PList_next(matches, &iterM))) {
          if (EN_id(match) != EN_id(vert)) {
            if (classType(match) == 0) {
              match_gEnt = classId(match);
            }
          }
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
#if (SIMMODSUITE_MAJOR_VERSION <= 14 && SIMMODSUITE_MINOR_VERSION < 190921) || \
    (SIMMODSUITE_MAJOR_VERSION >= 16 && SIMMODSUITE_MINOR_VERSION >= 210623)
      pPList matches = EN_getMatchingEnts(edge, NULL);
#else
      pPList matches = EN_getMatchingEnts(edge, NULL, 0);
#endif
      void *iterM = 0;
      pEdge match;
      if (matches) {
        assert(PList_size(matches) > 1);
        while ((match = (pEdge)PList_next(matches, &iterM))) {
          if (EN_id(match) != EN_id(edge)) {
            if (classType(match) == 1) {
              match_gEnt = classId(match);
            }
          }
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
#if (SIMMODSUITE_MAJOR_VERSION <= 14 && SIMMODSUITE_MINOR_VERSION < 190921) || \
    (SIMMODSUITE_MAJOR_VERSION >= 16 && SIMMODSUITE_MINOR_VERSION >= 210623)
      pPList matches = EN_getMatchingEnts(face, NULL);
#else
      pPList matches = EN_getMatchingEnts(face, NULL, 0);
#endif
      void *iterM = 0;
      pFace match;
      if (matches) {
        assert(PList_size(matches) > 1);
        while ((match = (pFace)PList_next(matches, &iterM))) {
          if (EN_id(match) != EN_id(face)) {
            if (classType(match) == 2) {
              match_gEnt = classId(match);
            }
          }
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
  EIter edges = M_edgeIter(m);
  pEdge edge;
  while ((edge = (pEdge) EIter_next(edges))) {
    for(int j=1; j>=0; --j) {
      vtx = E_vertex(edge, j);
      double xyz[3];
      V_coord(vtx,xyz);
      ent_nodes[1].push_back(EN_id(vtx));
    }
    ent_class_ids[1].push_back(classId(edge));

#if (SIMMODSUITE_MAJOR_VERSION <= 14 && SIMMODSUITE_MINOR_VERSION < 190921) || \
    (SIMMODSUITE_MAJOR_VERSION >= 16 && SIMMODSUITE_MINOR_VERSION >= 210623)
    pPList matches = EN_getMatchingEnts(edge, NULL);
#else
    pPList matches = EN_getMatchingEnts(edge, NULL, 0);
#endif
    void *iterM = 0;
    pEdge match;
    count_matches = 0;
    if (matches) {
      assert (PList_size(matches) > 1);
      while ((match = (pEdge)PList_next(matches, &iterM))) {
        if (EN_id(match) != EN_id(edge)) {
          ++count_matches;
          if (count_matches > 1) {
            Omega_h_fail("Error:matches per entity > 1\n");
          }

          leaf_idxs[1].push_back(EN_id(edge));
          root_idxs[1].push_back(EN_id(match));
          root_ranks[1].push_back(0);
        }
      }
      PList_delete(matches);
    }
  }
  EIter_delete(edges);
  const int numFaces = M_numFaces(m);
  ent_nodes[2].reserve(numFaces*3);
  ent_class_ids[2].reserve(numFaces);
  FIter faces = M_faceIter(m);
  pFace face;
  while ((face = (pFace) FIter_next(faces))) {
    pPList verts = F_vertices(face,1);
    assert(PList_size(verts) == 3);
    void *iter = 0; // must initialize to 0
    while ((vtx = (pVertex)PList_next(verts, &iter))) {
      ent_nodes[2].push_back(EN_id(vtx));
    }
    PList_delete(verts);
    ent_class_ids[2].push_back(classId(face));

#if (SIMMODSUITE_MAJOR_VERSION <= 14 && SIMMODSUITE_MINOR_VERSION < 190921) || \
    (SIMMODSUITE_MAJOR_VERSION >= 16 && SIMMODSUITE_MINOR_VERSION >= 210623)
    pPList matches = EN_getMatchingEnts(face, NULL);
#else
    pPList matches = EN_getMatchingEnts(face, NULL, 0);
#endif
    void *iterM = 0;
    pFace match;
    count_matches = 0;
    if (matches) {
      assert (PList_size(matches) > 1);
      while ((match = (pFace)PList_next(matches, &iterM))) {
        if ((PList_size(matches)>1) && (EN_id(match) != EN_id(face))) {

          ++count_matches;
          if (count_matches > 1) {
            Omega_h_fail("Error:matches per entity > 1\n");
          }

          leaf_idxs[2].push_back(EN_id(face));
          root_idxs[2].push_back(EN_id(match));
          root_ranks[2].push_back(0);
        }
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
#if (SIMMODSUITE_MAJOR_VERSION <= 14 && SIMMODSUITE_MINOR_VERSION < 190921) || \
    (SIMMODSUITE_MAJOR_VERSION >= 16 && SIMMODSUITE_MINOR_VERSION >= 210623)
    pPList matches = EN_getMatchingEnts(rgn, NULL);
#else
    pPList matches = EN_getMatchingEnts(rgn, NULL, 0);
#endif
    count_matches = 0;
    if (matches) {
      Omega_h_fail("region matches found\n");
    }
    PList_delete(matches);

    pPList verts = R_vertices(rgn, 1);
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
    mesh->add_tag(ent_dim, "global", 1, GOs(ndim_ents, 0, 1));
  }
  //
  if (!comm->reduce_and(is_sorted(vert_globals))) {
    reorder_by_globals(mesh);
  }
  mesh->add_coords(host_coords.write());
  //

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
  return;
}

void matchRead(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm, Mesh *mesh, I8 is_in) {
  mesh->set_matched(1);
  if (is_in) {
    MS_init();
    SimModel_start();
    pNativeModel nm = NULL;
    pProgress p = NULL;
    pGModel g = GM_load(mdl_fname.c_str(), nm, p);
    pMesh sm = M_load(mesh_fname.c_str(), g, p);

    mesh->set_comm(comm);
    meshsim::read_matchInternal(sm, mesh, g, comm);

    M_release(sm);
    GM_release(g);
    SimModel_stop();
    MS_exit();
    Sim_logOff();
  }
  return;
}

}  // namespace meshsim

}  // end namespace Omega_h
