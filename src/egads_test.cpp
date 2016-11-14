#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <egads.h>

#ifdef __clang__
#pragma clang diagnostic pop
#endif


#include <cassert>
#include <array>

#define CALL(f) assert(EGADS_SUCCESS == (f))

static char const* get_oclass_name(int oclass) {
  switch (oclass) {
    case NODE: return "node";
    case EDGE: return "edge";
    case LOOP: return "loop";
    case FACE: return "face";
    case SHELL: return "shell";
    case BODY: return "body";
    case MODEL: return "model";
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
  CALL(EG_getTopology(vert, &ref, &oclass, &mtype,
        data.data(), &nchild, &children, &senses));
  return data;
}

int main(int argc, char** argv) {
  assert(argc == 2);
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
  CALL(EG_getTopology(model, &model_geom, &model_oclass, &model_mtype,
      nullptr, &nbodies, &bodies, &body_senses));
  printf("model oclass \"%s\" mtype %d has %d bodies\n",
      get_oclass_name(model_oclass), model_mtype, nbodies);
  auto body = bodies[0];
  int nfaces;
  ego* faces;
  CALL(EG_getBodyTopos(body, nullptr, FACE, &nfaces, &faces));
  printf("first body has %d faces\n", nfaces);
  EG_free(faces);
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
    printf("point %d at %f %f %f\n", i, pt[0], pt[1], pt[2]);
  }
  EG_free(verts);
  CALL(EG_deleteObject(model));
  CALL(EG_close(context));
}
