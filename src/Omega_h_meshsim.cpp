#include "Omega_h_file.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_adj.hpp"

#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"
#include "SimDiscrete.h"

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

void read_internal(pMesh m, Mesh* mesh) {

  (void)mesh;

  RIter regions = M_regionIter(m);
  LO count_tet = 0;
  LO count_hex = 0;
  LO count_wedge = 0;
  LO count_pyramid = 0;
  pRegion rgn;
  while ((rgn = (pRegion) RIter_next(regions))) {
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

  FIter faces = M_faceIter(m);
  pFace face;
  int count_tri = 0;
  int count_quad = 0;
  while ((face = (pFace) FIter_next(faces))) {
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

  bool is_simplex = 0;
  bool is_hypercube = 0;
  if (count_hex == 0 && count_wedge == 0 && count_pyramid == 0) {
    if (count_tet == 0) {
      if (count_tri > 0) {
        is_simplex = 1;
      }
    }
    else if (count_tet > 0) {
      is_simplex = 1;
    }
    else {
      Omega_h_fail ("Invaild topology type\n");
    }
  }

  if (count_tet == 0 && count_wedge == 0 && count_pyramid == 0) {
    if (count_hex == 0) {
      if (count_quad > 0) {
        is_hypercube = 1;
      }
    }
    else if (count_hex > 0) {
      is_hypercube = 1;
    }
    else {
      Omega_h_fail ("Invaild topology type\n");
    }
  }

  fprintf(stderr, "tet=%d, hex=%d, wedge=%d, pyramid=%d\n",
         count_tet, count_hex, count_wedge, count_pyramid);
  fprintf(stderr, "tri=%d, quad=%d\n", count_tri, count_quad);

  const int numVtx = M_numVertices(m);
  const int numEdges = M_numEdges(m);
  const int numFaces = M_numFaces(m);
  const int numRegions = M_numRegions(m);

  std::vector<int> rgn_vertices[4];
  std::vector<int> face_vertices[2];
  std::vector<int> edge_vertices[1];
  std::vector<int> ent_class_ids[4];
  std::vector<int> ent_class_dim[4];

  ent_class_ids[0].reserve(numVtx);
  ent_class_dim[0].reserve(numVtx);
  ent_class_ids[1].reserve(numEdges);
  ent_class_dim[1].reserve(numEdges);
  ent_class_ids[2].reserve(numFaces);
  ent_class_dim[2].reserve(numFaces);
  ent_class_ids[3].reserve(numRegions);
  ent_class_dim[3].reserve(numRegions);

  Int max_dim;
  if (numRegions) {
    max_dim = 3;
  } else if (numFaces) {
    max_dim = 2;
  } else if (numEdges) {
    max_dim = 1;
  } else {
    Omega_h_fail("There were no Elements of dimension higher than zero!\n");
  }

  HostWrite<Real> host_coords(numVtx*max_dim);
  VIter vertices = M_vertexIter(m);
  pVertex vtx;
  LO v = 0;
  while ((vtx = (pVertex) VIter_next(vertices))) {
    double xyz[3];
    V_coord(vtx,xyz);
    if( max_dim < 3 && xyz[2] != 0 )
      Omega_h_fail("The z coordinate must be zero for a 2d mesh!\n");
    for(int j=0; j<max_dim; j++) {
      host_coords[v * max_dim + j] = xyz[j];
    }
    ent_class_ids[0].push_back(classId(vtx));
    ent_class_dim[0].push_back(classType(vtx));
    ++v;
  }
  VIter_delete(vertices);

  HostWrite<LO> host_class_ids_vtx(numVtx);
  HostWrite<I8> host_class_dim_vtx(numVtx);
  for (int i = 0; i < numVtx; ++i) {
    host_class_ids_vtx[i] = ent_class_ids[0][static_cast<std::size_t>(i)];
    host_class_dim_vtx[i] = ent_class_dim[0][static_cast<std::size_t>(i)];
  }

  mesh->set_dim(max_dim);
  if (is_simplex || is_hypercube) {
    if (is_simplex) {
      mesh->set_family(OMEGA_H_SIMPLEX);
    }
    else if (is_hypercube){
      mesh->set_family(OMEGA_H_HYPERCUBE);
    }
    mesh->set_verts(numVtx);
    mesh->add_coords(host_coords.write());
    mesh->add_tag<ClassId>(0, "class_id", 1,
                           Read<ClassId>(host_class_ids_vtx.write()));
    mesh->add_tag<I8>(0, "class_dim", 1,
                      Read<I8>(host_class_dim_vtx.write()));
  }
  else {
    mesh->set_family(OMEGA_H_MIXED);
    mesh->set_verts_type(numVtx);
    mesh->add_coords_mix(host_coords.write());
    mesh->add_tag<ClassId>(Topo_type::vertex, "class_id", 1,
                           Read<ClassId>(host_class_ids_vtx.write()));
    mesh->add_tag<I8>(Topo_type::vertex, "class_dim", 1,
                      Read<I8>(host_class_dim_vtx.write()));
  }

  edge_vertices[0].reserve(numEdges*2);
  EIter edges = M_edgeIter(m);
  pEdge edge;
  int count_edge = 0;
  while ((edge = (pEdge) EIter_next(edges))) {
    double xyz[3];
    count_edge += 1;
    for(int j=0; j<2; ++j) {
      vtx = E_vertex(edge,j);
      edge_vertices[0].push_back(EN_id(vtx));
      V_coord(vtx,xyz);
    }
    ent_class_ids[1].push_back(classId(edge));
    ent_class_dim[1].push_back(classType(edge));
  }
  EIter_delete(edges);
  
  HostWrite<LO> host_class_ids_edge(numEdges);
  HostWrite<I8> host_class_dim_edge(numEdges);
  for (int i = 0; i < numEdges; ++i) {
    host_class_ids_edge[i] = ent_class_ids[1][static_cast<std::size_t>(i)];
    host_class_dim_edge[i] = ent_class_dim[1][static_cast<std::size_t>(i)];
  }

  HostWrite<LO> host_e2v(numEdges*2);
  for (Int i = 0; i < numEdges; ++i) {
    for (Int j = 0; j < 2; ++j) {
      host_e2v[i*2 + j] =
          edge_vertices[0][static_cast<std::size_t>(i*2 + j)];
    }
  }
  auto ev2v = Read<LO>(host_e2v.write());
  if (is_simplex || is_hypercube) {
    mesh->set_ents(1, Adj(ev2v));
    mesh->add_tag<ClassId>(1, "class_id", 1,
                           Read<ClassId>(host_class_ids_edge.write()));
    mesh->add_tag<I8>(1, "class_dim", 1,
                      Read<I8>(host_class_dim_edge.write()));
  }
  else {
    mesh->set_ents(Topo_type::edge, Topo_type::vertex, Adj(ev2v));
    mesh->add_tag<ClassId>(Topo_type::edge, "class_id", 1,
                           Read<ClassId>(host_class_ids_edge.write()));
    mesh->add_tag<I8>(Topo_type::edge, "class_dim", 1,
                      Read<I8>(host_class_dim_edge.write()));
  }

  face_vertices[0].reserve(count_tri*3);
  face_vertices[1].reserve(count_quad*4);
  std::vector<int> face_class_ids[2];
  std::vector<int> face_class_dim[2];
  face_class_ids[0].reserve(count_tri);
  face_class_dim[0].reserve(count_tri);
  face_class_ids[1].reserve(count_quad);
  face_class_dim[1].reserve(count_quad);

  faces = M_faceIter(m);
  while ((face = (pFace) FIter_next(faces))) {
    if (F_numEdges(face) == 3) {
      pVertex tri_vertex;
      pPList tri_vertices = F_vertices(face,1);
      assert (PList_size(tri_vertices) == 3);
      void *iter = 0;
      while ((tri_vertex = (pVertex) PList_next(tri_vertices, &iter))) {
        face_vertices[0].push_back(EN_id(tri_vertex));
      }
      PList_delete(tri_vertices);
      face_class_ids[0].push_back(classId(face));
      face_class_dim[0].push_back(classType(face));
    }
    else if (F_numEdges(face) == 4) {
      pVertex quad_vertex;
      pPList quad_vertices = F_vertices(face,1);
      assert (PList_size(quad_vertices) == 4);
      void *iter = 0;
      while ((quad_vertex = (pVertex) PList_next(quad_vertices, &iter))) {
        face_vertices[1].push_back(EN_id(quad_vertex));
      }
      PList_delete(quad_vertices);
      face_class_ids[1].push_back(classId(face));
      face_class_dim[1].push_back(classType(face));
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
    ent_class_ids[2].push_back(classId(face));
    ent_class_dim[2].push_back(classType(face));
  }
  FIter_delete(faces);

  HostWrite<LO> host_class_ids_face(numFaces);
  HostWrite<I8> host_class_dim_face(numFaces);
  for (int i = 0; i < numFaces; ++i) {
    host_class_ids_face[i] = ent_class_ids[2][static_cast<std::size_t>(i)];
    host_class_dim_face[i] = ent_class_dim[2][static_cast<std::size_t>(i)];
  }
  HostWrite<LO> host_class_ids_tri(count_tri);
  HostWrite<I8> host_class_dim_tri(count_tri);
  for (int i = 0; i < count_tri; ++i) {
    host_class_ids_tri[i] = face_class_ids[0][static_cast<std::size_t>(i)];
    host_class_dim_tri[i] = face_class_dim[0][static_cast<std::size_t>(i)];
  }
  HostWrite<LO> host_class_ids_quad(count_quad);
  HostWrite<I8> host_class_dim_quad(count_quad);
  for (int i = 0; i < count_quad; ++i) {
    host_class_ids_quad[i] = face_class_ids[1][static_cast<std::size_t>(i)];
    host_class_dim_quad[i] = face_class_dim[1][static_cast<std::size_t>(i)];
  }

  Adj edge2vert;
  Adj vert2edge;
  if (is_simplex || is_hypercube) {
    edge2vert = mesh->get_adj(1, 0);
    vert2edge = mesh->ask_up(0, 1);
  }
  else {
    edge2vert = mesh->get_adj(Topo_type::edge, Topo_type::vertex);
    vert2edge = mesh->ask_up(Topo_type::vertex, Topo_type::edge);
  }
  HostWrite<LO> host_tri2verts(count_tri*3);
  for (Int i = 0; i < count_tri; ++i) {
    for (Int j = 0; j < 3; ++j) {
      host_tri2verts[i*3 + j] =
          face_vertices[0][static_cast<std::size_t>(i*3 + j)];
    }
  }
  auto tri2verts = Read<LO>(host_tri2verts.write());
  Adj down;
  if (is_simplex) {
    down = reflect_down(tri2verts, edge2vert.ab2b, vert2edge,
                        OMEGA_H_SIMPLEX, 2, 1);
    mesh->set_ents(2, down);
    mesh->add_tag<ClassId>(2, "class_id", 1,
                           Read<ClassId>(host_class_ids_face.write()));
    mesh->add_tag<I8>(2, "class_dim", 1,
                      Read<I8>(host_class_dim_face.write()));
  }
  else if (is_hypercube) {
    //empty since a quad/hex mesh has no triangles and to avoid dropping into
    //the 'else'
  }
  else {
    down = reflect_down(tri2verts, edge2vert.ab2b, vert2edge,
                        Topo_type::triangle, Topo_type::edge);
    mesh->set_ents(Topo_type::triangle, Topo_type::edge, down);
    mesh->add_tag<ClassId>(Topo_type::triangle, "class_id", 1,
                           Read<ClassId>(host_class_ids_tri.write()));
    mesh->add_tag<I8>(Topo_type::triangle, "class_dim", 1,
                      Read<I8>(host_class_dim_tri.write()));
  }

  HostWrite<LO> host_quad2verts(count_quad*4);
  for (Int i = 0; i < count_quad; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_quad2verts[i*4 + j] =
          face_vertices[1][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto quad2verts = Read<LO>(host_quad2verts.write());

  if (is_hypercube) {
    down = reflect_down(quad2verts, edge2vert.ab2b, vert2edge,
                        OMEGA_H_HYPERCUBE, 2, 1);
    mesh->set_ents(2, down);
    mesh->add_tag<ClassId>(2, "class_id", 1,
                           Read<ClassId>(host_class_ids_face.write()));
    mesh->add_tag<I8>(2, "class_dim", 1,
                      Read<I8>(host_class_dim_face.write()));
  }
  else if (is_simplex) {
    //empty to avoid dropping into the 'else'
  }
  else {
    down = reflect_down(quad2verts, edge2vert.ab2b, vert2edge,
                      Topo_type::quadrilateral, Topo_type::edge);
    mesh->set_ents(Topo_type::quadrilateral, Topo_type::edge, down);
    mesh->add_tag<ClassId>(Topo_type::quadrilateral, "class_id", 1,
                           Read<ClassId>(host_class_ids_quad.write()));
    mesh->add_tag<I8>(Topo_type::quadrilateral, "class_dim", 1,
                      Read<I8>(host_class_dim_quad.write()));
  }

  if (!(count_tet == 0 && count_hex == 0 && count_wedge == 0 && 
        count_pyramid == 0)) {
    rgn_vertices[0].reserve(count_tet*4);
    rgn_vertices[1].reserve(count_hex*8);
    rgn_vertices[2].reserve(count_wedge*6);
    rgn_vertices[3].reserve(count_pyramid*5);
    std::vector<int> rgn_class_ids[4];
    std::vector<int> rgn_class_dim[4];
    rgn_class_ids[0].reserve(count_tet);
    rgn_class_dim[0].reserve(count_tet);
    rgn_class_ids[1].reserve(count_hex);
    rgn_class_dim[1].reserve(count_hex);
    rgn_class_ids[2].reserve(count_wedge);
    rgn_class_dim[2].reserve(count_wedge);
    rgn_class_ids[3].reserve(count_pyramid);
    rgn_class_dim[3].reserve(count_pyramid);

    regions = M_regionIter(m);
    while ((rgn = (pRegion) RIter_next(regions))) {
      if (R_topoType(rgn) == Rtet) {
        pVertex vert;
        pPList verts = R_vertices(rgn,1);
        assert (PList_size(verts) == 4);
        void *iter = 0;
        while ((vert = (pVertex) PList_next(verts, &iter))) {
          rgn_vertices[0].push_back(EN_id(vert));
        }
        PList_delete(verts);
        rgn_class_ids[0].push_back(classId(rgn));
        rgn_class_dim[0].push_back(classType(rgn));
      }
      else if (R_topoType(rgn) == Rhex) {
        pVertex vert;
        pPList verts = R_vertices(rgn,1);
        assert (PList_size(verts) == 8);
        void *iter = 0;
        while ((vert = (pVertex) PList_next(verts, &iter))) {
          rgn_vertices[1].push_back(EN_id(vert));
        }
        PList_delete(verts);
        rgn_class_ids[1].push_back(classId(rgn));
        rgn_class_dim[1].push_back(classType(rgn));
      }
      else if (R_topoType(rgn) == Rwedge) {
        pVertex vert;
        pPList verts = R_vertices(rgn,1);
        assert (PList_size(verts) == 6);
        void *iter = 0;
        while ((vert = (pVertex) PList_next(verts, &iter))) {
          rgn_vertices[2].push_back(EN_id(vert));
        }
        PList_delete(verts);
        rgn_class_ids[2].push_back(classId(rgn));
        rgn_class_dim[2].push_back(classType(rgn));
      }
      else if (R_topoType(rgn) == Rpyramid) {
        pVertex vert;
        pPList verts = R_vertices(rgn,1);
        assert (PList_size(verts) == 5);
        void *iter = 0;
        while ((vert = (pVertex) PList_next(verts, &iter))) {
          rgn_vertices[3].push_back(EN_id(vert));
        }
        PList_delete(verts);
        rgn_class_ids[3].push_back(classId(rgn));
        rgn_class_dim[3].push_back(classType(rgn));
      }
      else {
        Omega_h_fail ("Region is not tet, hex, wedge, or pyramid \n");
      }
      ent_class_ids[3].push_back(classId(rgn));
      ent_class_dim[3].push_back(classType(rgn));
    }
    RIter_delete(regions);

    HostWrite<LO> host_class_ids_rgn(numRegions);
    HostWrite<I8> host_class_dim_rgn(numRegions);
    for (int i = 0; i < numRegions; ++i) {
      host_class_ids_rgn[i] = ent_class_ids[3][static_cast<std::size_t>(i)];
      host_class_dim_rgn[i] = ent_class_dim[3][static_cast<std::size_t>(i)];
    }
    HostWrite<LO> host_class_ids_tet(count_tet);
    HostWrite<I8> host_class_dim_tet(count_tet);
    for (int i = 0; i < count_tet; ++i) {
      host_class_ids_tet[i] = rgn_class_ids[0][static_cast<std::size_t>(i)];
      host_class_dim_tet[i] = rgn_class_dim[0][static_cast<std::size_t>(i)];
    }
    HostWrite<LO> host_class_ids_hex(count_hex);
    HostWrite<I8> host_class_dim_hex(count_hex);
    for (int i = 0; i < count_hex; ++i) {
      host_class_ids_hex[i] = rgn_class_ids[1][static_cast<std::size_t>(i)];
      host_class_dim_hex[i] = rgn_class_dim[1][static_cast<std::size_t>(i)];
    }
    HostWrite<LO> host_class_ids_wedge(count_wedge);
    HostWrite<I8> host_class_dim_wedge(count_wedge);
    for (int i = 0; i < count_wedge; ++i) {
      host_class_ids_wedge[i] = rgn_class_ids[2][static_cast<std::size_t>(i)];
      host_class_dim_wedge[i] = rgn_class_dim[2][static_cast<std::size_t>(i)];
    }
    HostWrite<LO> host_class_ids_pyramid(count_pyramid);
    HostWrite<I8> host_class_dim_pyramid(count_pyramid);
    for (int i = 0; i < count_pyramid; ++i) {
      host_class_ids_pyramid[i] = rgn_class_ids[3][static_cast<std::size_t>(i)];
      host_class_dim_pyramid[i] = rgn_class_dim[3][static_cast<std::size_t>(i)];
    }

    Adj tri2vert;
    Adj vert2tri;
    Adj quad2vert;
    Adj vert2quad;
    if (is_simplex) {
      tri2vert = mesh->ask_down(2, 0);
      vert2tri = mesh->ask_up(0, 2);
    }
    else if (is_hypercube) {
      quad2vert = mesh->ask_down(2, 0);
      vert2quad = mesh->ask_up(0, 2);
    }
    else {
      tri2vert = mesh->ask_down(Topo_type::triangle, Topo_type::vertex);
      vert2tri = mesh->ask_up(Topo_type::vertex, Topo_type::triangle);
      quad2vert = mesh->ask_down(Topo_type::quadrilateral, Topo_type::vertex);
      vert2quad = mesh->ask_up(Topo_type::vertex, Topo_type::quadrilateral);
    }

    HostWrite<LO> host_tet2verts(count_tet*4);
    for (Int i = 0; i < count_tet; ++i) {
      for (Int j = 0; j < 4; ++j) {
        host_tet2verts[i*4 + j] =
          rgn_vertices[0][static_cast<std::size_t>(i*4 + j)];
      }
    }
    auto tet2verts = Read<LO>(host_tet2verts.write());
    if (is_simplex) {
      down = reflect_down(tet2verts, tri2vert.ab2b, vert2tri,
          OMEGA_H_SIMPLEX, 3, 2);
      mesh->set_ents(3, down);
      mesh->add_tag<ClassId>(3, "class_id", 1,
          Read<ClassId>(host_class_ids_rgn.write()));
      mesh->add_tag<I8>(3, "class_dim", 1,
          Read<I8>(host_class_dim_rgn.write()));
    }
    else if (is_hypercube) {
      //empty to avoid dropping into the 'else'
    }
    else {
      down = reflect_down(tet2verts, tri2vert.ab2b, vert2tri,
          Topo_type::tetrahedron, Topo_type::triangle);
      mesh->set_ents(Topo_type::tetrahedron, Topo_type::triangle, down);
      mesh->add_tag<ClassId>(Topo_type::tetrahedron, "class_id", 1,
          Read<ClassId>(host_class_ids_tet.write()));
      mesh->add_tag<I8>(Topo_type::tetrahedron, "class_dim", 1,
          Read<I8>(host_class_dim_tet.write()));
    }

    HostWrite<LO> host_hex2verts(count_hex*8);
    for (Int i = 0; i < count_hex; ++i) {
      for (Int j = 0; j < 8; ++j) {
        host_hex2verts[i*8 + j] =
          rgn_vertices[1][static_cast<std::size_t>(i*8 + j)];
      }
    }
    auto hex2verts = Read<LO>(host_hex2verts.write());

    if (is_hypercube) {
      down = reflect_down(hex2verts, quad2vert.ab2b, vert2quad,
          OMEGA_H_HYPERCUBE, 3, 2);
      mesh->set_ents(3, down);
      mesh->add_tag<ClassId>(3, "class_id", 1,
          Read<ClassId>(host_class_ids_rgn.write()));
      mesh->add_tag<I8>(3, "class_dim", 1,
          Read<I8>(host_class_dim_rgn.write()));
    }
    else if (is_simplex) {
      //empty to avoid dropping into the 'else'
    }
    else {
      down = reflect_down(hex2verts, quad2vert.ab2b, vert2quad,
          Topo_type::hexahedron, Topo_type::quadrilateral);
      mesh->set_ents(Topo_type::hexahedron, Topo_type::quadrilateral, down);
      mesh->add_tag<ClassId>(Topo_type::hexahedron, "class_id", 1,
          Read<ClassId>(host_class_ids_hex.write()));
      mesh->add_tag<I8>(Topo_type::hexahedron, "class_dim", 1,
          Read<I8>(host_class_dim_hex.write()));
    }

    HostWrite<LO> host_wedge2verts(count_wedge*6);
    for (Int i = 0; i < count_wedge; ++i) {
      for (Int j = 0; j < 6; ++j) {
        host_wedge2verts[i*6 + j] =
          rgn_vertices[2][static_cast<std::size_t>(i*6 + j)];
      }
    }
    auto wedge2verts = Read<LO>(host_wedge2verts.write());
    down = reflect_down(wedge2verts, quad2vert.ab2b, vert2quad,
        Topo_type::wedge, Topo_type::quadrilateral);
    if ((!is_simplex) && (!is_hypercube)) {
      mesh->set_ents(Topo_type::wedge, Topo_type::quadrilateral, down);
    }

    down = reflect_down(wedge2verts, tri2vert.ab2b, vert2tri,
        Topo_type::wedge, Topo_type::triangle);
    if ((!is_simplex) && (!is_hypercube)) {
      mesh->set_ents(Topo_type::wedge, Topo_type::triangle, down);
      mesh->add_tag<ClassId>(Topo_type::wedge, "class_id", 1,
          Read<ClassId>(host_class_ids_wedge.write()));
      mesh->add_tag<I8>(Topo_type::wedge, "class_dim", 1,
          Read<I8>(host_class_dim_wedge.write()));
    }

    HostWrite<LO> host_pyramid2verts(count_pyramid*5);
    for (Int i = 0; i < count_pyramid; ++i) {
      for (Int j = 0; j < 5; ++j) {
        host_pyramid2verts[i*5 + j] =
          rgn_vertices[3][static_cast<std::size_t>(i*5 + j)];
      }
    }
    auto pyramid2verts = Read<LO>(host_pyramid2verts.write());
    down = reflect_down(pyramid2verts, tri2vert.ab2b, vert2tri,
        Topo_type::pyramid, Topo_type::triangle);
    if ((!is_simplex) && (!is_hypercube)) {
      mesh->set_ents(Topo_type::pyramid, Topo_type::triangle, down);
    }

    down = reflect_down(pyramid2verts, quad2vert.ab2b, vert2quad,
        Topo_type::pyramid, Topo_type::quadrilateral);
    if ((!is_simplex) && (!is_hypercube)) {
      mesh->set_ents(Topo_type::pyramid, Topo_type::quadrilateral, down);
      mesh->add_tag<ClassId>(Topo_type::pyramid, "class_id", 1,
          Read<ClassId>(host_class_ids_pyramid.write()));
      mesh->add_tag<I8>(Topo_type::pyramid, "class_dim", 1,
          Read<I8>(host_class_dim_pyramid.write()));
    }
  }

  return;
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
  mesh.set_comm(comm);
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  meshsim::read_internal(m, &mesh);
  M_release(m);
  GM_release(g);
  SimDiscrete_stop(0);
  SimModel_stop();
  return mesh;
}

}  // namespace meshsim

}  // end namespace Omega_h
