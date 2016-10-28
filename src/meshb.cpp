#include <Omega_h.hpp>
#include <libmeshb7.h>

namespace Omega_h {

namespace meshb {

template <int version>
struct VersionTypes;

template <>
struct VersionTypes<1> {
  using Index = std::int32_t;
  using Real = float;
  using Line = std::int32_t;
};

template <>
struct VersionTypes<2> {
  using Index = std::int32_t;
  using Real = double;
  using Line = std::int32_t;
};

template <>
struct VersionTypes<3> {
  using Index = std::int32_t;
  using Real = double;
  using Line = std::int64_t;
};

template <>
struct VersionTypes<4> {
  using Index = std::int64_t;
  using Real = double;
  using Line = std::int64_t;
};

static void safe_goto(long long file, GmfKwdCod key) {
  auto ret = GmfGotoKwd(file, key);
  OMEGA_H_CHECK(ret);
}

template <int version>
static void read_meshb_version(
    Mesh* mesh, long long file, int ver, int dim) {
  using GmfIndex = typename VersionTypes<version>::Index;
  using GmfReal = typename VersionTypes<version>::Real;
  using GmfLine = typename VersionTypes<version>::Line;
  OMEGA_H_CHECK(ver == version);
  GmfIndex nverts = GmfStatKwd(file, GmfVertices);
  static GmfKwdCod const simplex_kwds[4] = {
      GmfVertices, GmfEdges, GmfTriangles, GmfTetrahedra};
  safe_goto(file, GmfVertices);
  auto coords = HostWrite<Real>(nverts * dim);
  for (GmfLine i = 0; i < nverts; ++i) {
    GmfIndex ref;
    Few<GmfReal, 3> tmp;
    if (dim == 2) GmfGetLin(file, GmfVertices, &tmp[0], &tmp[1], &ref);
    if (dim == 3) GmfGetLin(file, GmfVertices, &tmp[0], &tmp[1], &tmp[2], &ref);
    OMEGA_H_CHECK(ref == i + 1);
    for (int j = 0; j < dim; ++j) coords[i * dim + j] = tmp[j];
  }
  auto side_kwd = simplex_kwds[dim - 1];
  GmfIndex nsides = GmfStatKwd(file, side_kwd);
  auto sides2verts = HostWrite<LO>(nsides * dim);
  auto sides2class_id = HostWrite<I32>(nsides);
  if (nsides) {
    safe_goto(file, side_kwd);
    for (GmfLine i = 0; i < nsides; ++i) {
      GmfIndex class_id;
      Few<GmfIndex, 3> tmp;
      if (dim == 2) GmfGetLin(file, side_kwd, &tmp[0], &tmp[1], &class_id);
      if (dim == 3)
        GmfGetLin(file, side_kwd, &tmp[0], &tmp[1], &tmp[2], &class_id);
      for (int j = 0; j < dim; ++j) sides2verts[i * dim + j] = tmp[j] - 1;
      sides2class_id[i] = class_id;
    }
  }
  auto elem_kwd = simplex_kwds[dim];
  GmfIndex nelems = GmfStatKwd(file, elem_kwd);
  safe_goto(file, elem_kwd);
  auto elems2verts = HostWrite<LO>(nelems * (dim + 1));
  for (GmfLine i = 0; i < nelems; ++i) {
    GmfIndex ref;
    Few<GmfIndex, 4> tmp;
    if (dim == 2) GmfGetLin(file, elem_kwd, &tmp[0], &tmp[1], &tmp[2], &ref);
    if (dim == 3)
      GmfGetLin(file, elem_kwd, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &ref);
    OMEGA_H_CHECK(ref == i + 1);
    for (int j = 0; j <= dim; ++j) elems2verts[i * (dim + 1) + j] = tmp[j] - 1;
  }
  build_from_elems2verts(mesh, dim, LOs(elems2verts.write()), nverts);
  mesh->add_tag(VERT, "coordinates", dim, OMEGA_H_LINEAR_INTERP,
      OMEGA_H_DO_OUTPUT, Reals(coords.write()));
  GmfCloseMesh(file);
}

void read(Mesh* mesh, const char* filepath) {
  int ver, dim;
  auto file = GmfOpenMesh(filepath, GmfRead, &ver, &dim);
  if (!file) {
    Omega_h_fail("could not open Meshb file %s\n", filepath);
  }
  OMEGA_H_CHECK(dim == 2 || dim == 3);
  switch (ver) {
    case 1:
      read_meshb_version<1>(mesh, file, ver, dim);
      return;
    case 2:
      read_meshb_version<2>(mesh, file, ver, dim);
      return;
    case 3:
      read_meshb_version<3>(mesh, file, ver, dim);
      return;
    case 4:
      read_meshb_version<4>(mesh, file, ver, dim);
      return;
  }
  Omega_h_fail("unknown libMeshb version %d\n", ver);
}

}

}
