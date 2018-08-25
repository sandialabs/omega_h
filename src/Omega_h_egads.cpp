#include "Omega_h_egads.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_timer.hpp"

#include <Omega_h_fail.hpp>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <egads.h>

/* EGADS pollutes the global namespace with
 * macros like EDGE and FACE.
 * This is really bad.
 * I found a workaround for integer-valued macros,
 * but I hope EGADS stops doing this.
 */

enum EgadsObjectClass {
  EGADS_CONTXT = CONTXT,
  EGADS_TRANSFORM = TRANSFORM,
  EGADS_TESSELATION = TESSELLATION,
  EGADS_NIL = NIL,
  /*EGADS_EMPTY = EMPTY, not doing this one
   * because an EGADS error exists by the same name
   */
  EGADS_REFERENCE = REFERENCE,
  EGADS_PCURVE = PCURVE,
  EGADS_CURVE = CURVE,
  EGADS_SURFACE = SURFACE,
  EGADS_NODE = NODE,
  EGADS_EDGE = EDGE,
  EGADS_LOOP = LOOP,
  EGADS_FACE = FACE,
  EGADS_SHELL = SHELL,
  EGADS_BODY = BODY,
  EGADS_MODEL = MODEL
};

#undef CONTXT
#undef TRANSFORM
#undef TESSELLATION
#undef NIL
#undef EMPTY
#undef REFERENCE
#undef PCURVE
#undef CURVE
#undef SURFACE
#undef NODE
#undef EDGE
#undef LOOP
#undef FACE
#undef SHELL
#undef BODY
#undef MODEL

#ifdef __clang__
#pragma clang diagnostic pop
#endif

namespace Omega_h {

static void call_egads(
    int result, char const* code, char const* file, int line) {
  if (EGADS_SUCCESS == result) return;
  Omega_h_fail(
      "EGADS call %s returned %d at %s +%d\n", code, result, file, line);
}

#define CALL(f) call_egads((f), #f, __FILE__, __LINE__)

static int const dims2oclass[4] = {
    EGADS_NODE, EGADS_EDGE, EGADS_FACE, EGADS_BODY};

struct Egads {
  ego context;
  ego model;
  ego body;
  int counts[3];
  ego* entities[3];
  std::map<std::set<ego>, ego> classifier;
};

Egads* egads_load(std::string const& filename) {
  auto eg = new Egads;
  CALL(EG_open(&eg->context));
  CALL(EG_loadModel(eg->context, 0, filename.c_str(), &eg->model));
  ego model_geom;
  int model_oclass;
  int model_mtype;
  int nbodies;
  ego* bodies;
  int* body_senses;
  CALL(EG_getTopology(eg->model, &model_geom, &model_oclass, &model_mtype,
      nullptr, &nbodies, &bodies, &body_senses));
  OMEGA_H_CHECK(nbodies == 1);
  eg->body = bodies[0];
  for (int i = 0; i < 3; ++i) {
    CALL(EG_getBodyTopos(
        eg->body, nullptr, dims2oclass[i], &eg->counts[i], &eg->entities[i]));
  }
  // preprocess edge and vertex adjacency to faces
  for (int i = 0; i < 2; ++i) {
    std::vector<std::set<ego>> idxs2adj_faces(eg->counts[i]);
    for (int j = 0; j < eg->counts[2]; ++j) {
      auto face = eg->entities[2][j];
      int nadj_ents;
      ego* adj_ents;
      CALL(EG_getBodyTopos(
          eg->body, face, dims2oclass[i], &nadj_ents, &adj_ents));
      for (int k = 0; k < nadj_ents; ++k) {
        auto adj_ent = adj_ents[k];
        auto idx = EG_indexBodyTopo(eg->body, adj_ent) - 1;
        idxs2adj_faces[idx].insert(face);
      }
    }
    for (int j = 0; j < eg->counts[i]; ++j) {
      auto adj_faces = idxs2adj_faces[j];
      // HACK!: we have a really insane CAD model with nonsensical topology.
      // this essentially manifests as edges that are adjacent to only one
      // model face.
      // we actually want to just ignore these edges, so we won't create
      // classifier entries for them.
      if (adj_faces.size() == 1) continue;
      eg->classifier[adj_faces] = eg->entities[i][j];
    }
  }
  return eg;
}

static int get_dim(ego e) {
  ego ref;
  int oclass;
  int mtype;
  int nchild;
  ego* children;
  int* senses;
  CALL(EG_getTopology(
      e, &ref, &oclass, &mtype, nullptr, &nchild, &children, &senses));
  for (int i = 0; i <= 3; ++i)
    if (dims2oclass[i] == oclass) return i;
  return -1;
}

void egads_classify(Egads* eg, int nadj_faces, int const adj_face_ids[],
    int* class_dim, int* class_id) {
  std::set<ego> uniq_adj_faces;
  for (int i = 0; i < nadj_faces; ++i) {
    auto adj_face = eg->entities[2][adj_face_ids[i] - 1];
    uniq_adj_faces.insert(adj_face);
  }
  auto it = eg->classifier.find(uniq_adj_faces);
  if (it != eg->classifier.end()) {
    auto ent = it->second;
    *class_dim = get_dim(ent);
    *class_id = EG_indexBodyTopo(eg->body, ent);
  }
}

void egads_free(Egads* eg) {
  for (int i = 0; i < 3; ++i) {
    EG_free(eg->entities[i]);
  }
  CALL(EG_deleteObject(eg->model));
  CALL(EG_close(eg->context));
  delete eg;
}

void egads_reclassify(Mesh* mesh, Egads* eg) {
  OMEGA_H_CHECK(mesh->dim() == 3);
  auto face_class_dims = mesh->get_array<I8>(FACE, "class_dim");
  auto face_class_ids = mesh->get_array<ClassId>(FACE, "class_id");
  for (Int dim = 0; dim < 2; ++dim) {
    auto ents2faces = mesh->ask_up(dim, FACE);
    auto adj_class_dims = read(unmap(ents2faces.ab2b, face_class_dims, 1));
    auto keep_edges = each_eq_to(adj_class_dims, I8(2));
    auto ents2eq_faces = filter_graph(ents2faces, keep_edges);
    auto adj_eq_face_ids = unmap(ents2eq_faces.ab2b, face_class_ids, 1);
    auto host_a2ab = HostRead<LO>(ents2eq_faces.a2ab);
    auto host_face_ids = HostRead<LO>(adj_eq_face_ids);
    auto class_dims = mesh->get_array<I8>(dim, "class_dim");
    auto class_ids = mesh->get_array<ClassId>(dim, "class_id");
    auto host_class_dims = HostWrite<I8>(deep_copy(class_dims));
    auto host_class_ids = HostWrite<LO>(deep_copy(class_ids));
    for (LO i = 0; i < mesh->nents(dim); ++i) {
      auto b = host_a2ab[i];
      auto e = host_a2ab[i + 1];
      Int class_dim = host_class_dims[i];
      LO class_id = host_class_ids[i];
      egads_classify(
          eg, e - b, host_face_ids.data() + b, &class_dim, &class_id);
      host_class_dims[i] = I8(class_dim);
      host_class_ids[i] = class_id;
    }
    class_dims = Read<I8>(host_class_dims.write());
    class_ids = Read<LO>(host_class_ids.write());
    mesh->set_tag(dim, "class_id", class_ids);
    mesh->set_tag(dim, "class_dim", class_dims);
  }
}

static Vector<3> get_closest_point(ego g, Vector<3> in) {
  double ignored[2] = {};
  Vector<3> out = in;
  CALL(EG_invEvaluate(g, in.data(), ignored, out.data()));
  return out;
}

Reals egads_get_snap_warp(Mesh* mesh, Egads* eg, bool verbose) {
  OMEGA_H_CHECK(mesh->dim() == 3);
  if (verbose) std::cout << "Querying closest points for surface vertices...\n";
  auto t0 = now();
  auto class_dims = mesh->get_array<I8>(VERT, "class_dim");
  auto class_ids = mesh->get_array<ClassId>(VERT, "class_id");
  auto coords = mesh->coords();
  auto host_class_dims = HostRead<I8>(class_dims);
  auto host_class_ids = HostRead<LO>(class_ids);
  auto host_coords = HostRead<Real>(coords);
  auto host_warp = HostWrite<Real>(mesh->nverts() * 3);
  for (LO i = 0; i < mesh->nverts(); ++i) {
    auto a = get_vector<3>(host_coords, i);
    Int class_dim = host_class_dims[i];
    OMEGA_H_CHECK(class_dim >= 0);
    OMEGA_H_CHECK(class_dim <= 3);
    auto d = vector_3(0, 0, 0);
    if (0 < class_dim && class_dim < 3) {
      auto index = host_class_ids[i] - 1;
      OMEGA_H_CHECK(index >= 0);
      OMEGA_H_CHECK(index < eg->counts[class_dim]);
      auto g = eg->entities[class_dim][index];
      auto index2 = EG_indexBodyTopo(eg->body, g);
      OMEGA_H_CHECK(index2 == index + 1);
      auto b = get_closest_point(g, a);
      d = b - a;
    }
    set_vector(host_warp, i, d);
  }
  auto warp = Reals(host_warp.write());
  auto t1 = now();
  if (verbose) {
    std::cout << "Querying closest points for surface vertices took "
              << (t1 - t0) << " seconds\n";
  }
  return warp;
}

}  // namespace Omega_h
