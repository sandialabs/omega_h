#include "Omega_h_egads.hpp"
#include "internal.hpp"

#include <cassert>
#include <vector>
#include <map>
#include <set>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <egads.h>

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#define CALL(f) assert(EGADS_SUCCESS == (f))

namespace Omega_h {

static int const dims2oclass[4] = {
  NODE,
  EDGE,
  FACE,
  BODY
};

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
  CALL(EG_getTopology(eg->model, &model_geom, &model_oclass, &model_mtype, nullptr,
      &nbodies, &bodies, &body_senses));
  assert(nbodies == 1);
  eg->body = bodies[0];
  for (int i = 0; i < 3; ++i) {
    CALL(EG_getBodyTopos(eg->body, nullptr, dims2oclass[i],
         &eg->counts[i], &eg->entities[i]));
  }
  for (int i = 0; i < 2; ++i) {
    std::vector<std::set<ego>> idxs2adj_faces(eg->counts[i]);
    for (int j = 0; j < eg->counts[2]; ++j) {
      auto face = eg->entities[2][j];
      int nadj_ents;
      ego* adj_ents;
      CALL(EG_getBodyTopos(eg->body, face, dims2oclass[i],
            &nadj_ents, &adj_ents));
      for (int k = 0; k < nadj_ents; ++k) {
        auto adj_ent = adj_ents[k];
        auto idx = EG_indexBodyTopo(eg->body, adj_ent) - 1;
        idxs2adj_faces[idx].insert(face);
      }
    }
    for (int j = 0; j < eg->counts[i]; ++j) {
      auto adj_faces = idxs2adj_faces[j];
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
    if (dims2oclass[i] == oclass)
      return i;
  return -1;
}

void egads_classify(Egads* eg, int nadj_faces, int const adj_face_ids[],
    int* class_dim, int* class_id) {
  std::set<ego> uniq_adj_faces;
  for (int i = 0; i < nadj_faces; ++i) {
    auto adj_face = eg->entities[2][adj_face_ids[i] - 1];
    uniq_adj_faces.insert(adj_face);
  }
  auto ent = eg->classifier[uniq_adj_faces];
  *class_dim = get_dim(ent);
  *class_id = EG_indexBodyTopo(eg->body, ent);
}

void egads_free(Egads* eg) {
  for (int i = 0; i < 3; ++i) {
    EG_free(eg->entities[i]);
  }
  delete eg;
}

}
