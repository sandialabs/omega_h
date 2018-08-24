#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <egads.h>

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include <array>
#include <cassert>
#include <cmath>

#define CALL(f) OMEGA_H_CHECK(EGADS_SUCCESS == (f))

static char const* get_oclass_name(int oclass) {
  switch (oclass) {
    case NODE:
      return "node";
    case EDGE:
      return "edge";
    case LOOP:
      return "loop";
    case FACE:
      return "face";
    case SHELL:
      return "shell";
    case BODY:
      return "body";
    case MODEL:
      return "model";
  };
  return nullptr;
}

static std::array<double, 3> get_point(ego vert) {
  ego ref;
  int oclass;
  int mtype;
  std::array<double, 3> data;
  int nchild;
  ego* children;
  int* senses;
  CALL(EG_getTopology(
      vert, &ref, &oclass, &mtype, data.data(), &nchild, &children, &senses));
  return data;
}

static std::array<double, 3> get_closest_point(
    ego face, std::array<double, 3> in) {
  double ignored[2];
  std::array<double, 3> out = {{0.0, 0.0, 0.0}};
  auto ret = EG_invEvaluateGuess(face, in.data(), ignored, out.data());
  if (ret != EGADS_SUCCESS) {
    printf("EG_invEvaluateGuess failed with code %d\n", ret);
    ret = EG_invEvaluate(face, in.data(), ignored, out.data());
    if (ret != EGADS_SUCCESS) {
      printf("EG_invEvaluate also failed, with code %d\n", ret);
    }
  }
  return out;
}

static void print_closest_point(ego* faces, int i, std::array<double, 3> pt) {
  auto face = faces[i - 1];
  auto closest = get_closest_point(face, pt);
  printf("closest point to %f %f %f on face %d is %f %f %f\n", pt[0], pt[1],
      pt[2], i, closest[0], closest[1], closest[2]);
}

int main(int argc, char** argv) {
  OMEGA_H_CHECK(argc == 2);
  ego context;
  CALL(EG_open(&context));
  ego model;
  CALL(EG_loadModel(context, 0, argv[1], &model));
  ego model_geom;
  int model_oclass;
  int model_mtype;
  int nbodies;
  ego* bodies;
  int* body_senses;
  CALL(EG_getTopology(model, &model_geom, &model_oclass, &model_mtype, nullptr,
      &nbodies, &bodies, &body_senses));
  printf("model oclass \"%s\" mtype %d has %d bodies\n",
      get_oclass_name(model_oclass), model_mtype, nbodies);
  auto body = bodies[0];
  int nfaces;
  ego* faces;
  CALL(EG_getBodyTopos(body, nullptr, FACE, &nfaces, &faces));
  printf("first body has %d faces\n", nfaces);
  int nedges;
  ego* edges;
  CALL(EG_getBodyTopos(body, nullptr, EDGE, &nedges, &edges));
  printf("first body has %d edges\n", nedges);
  EG_free(edges);
  int nverts;
  ego* verts;
  CALL(EG_getBodyTopos(body, nullptr, NODE, &nverts, &verts));
  printf("first body has %d verts\n", nverts);
  for (int i = 0; i < nverts; ++i) {
    auto pt = get_point(verts[i]);
    printf("point %d at %f %f %f\n", i + 1, pt[0], pt[1], pt[2]);
  }
  EG_free(verts);
  print_closest_point(faces, 1, {{0.0, 0.75, 0.5}});
  print_closest_point(faces, 2, {{0.5, 0.5, 1.0}});
  print_closest_point(faces, 3, {{0.5, 1.0, 0.5}});
  print_closest_point(faces, 4, {{0.5, 0.5, 0.0}});
  print_closest_point(
      faces, 5, {{0.5 / std::sqrt(2), 0.5 / std::sqrt(2), 0.5}});
  print_closest_point(faces, 6, {{1.0, 0.5, 0.5}});
  print_closest_point(faces, 7, {{0.75, 0.0, 0.5}});
  EG_free(faces);
  CALL(EG_deleteObject(model));
  CALL(EG_close(context));
}
