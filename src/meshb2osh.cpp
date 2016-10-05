#include <libmeshb7.h>
#include <Omega_h.hpp>

#include <iostream>

namespace osh = Omega_h;

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

template <int version>
static void read_meshb_version(osh::Mesh* m, long long file, int ver, int dim) {
  using Index = typename VersionTypes<version>::Index;
  using Real = typename VersionTypes<version>::Real;
  using Line = typename VersionTypes<version>::Line;
  OMEGA_H_CHECK(ver == version);
  Index nverts = GmfStatKwd(file, GmfVertices);
  auto coords = osh::HostWrite<osh::Real>(nverts * dim);
  static GmfKwdCod const elem_kwds[4] = {
    GmfVertices,
    GmfEdges,
    GmfTriangles,
    GmfTetrahedra
  };
  auto elem_kwd = elem_kwds[dim];
  Index nelems = GmfStatKwd(file, elem_kwd);
  GmfGotoKwd(file, GmfVertices);
  for (Line i = 0; i < nverts; ++i) {
    Index ref;
    /* Omega_h doesn't instantiate its templates for `float`,
     * so `tmp` is here to temporarily store floats from version
     * 1 so they can be copied to doubles in Omega_h
     */
    osh::Few<Real, 3> tmp;
    if (dim == 2) GmfGetLin(file, GmfVertices, &tmp[0], &tmp[1], &ref);
    if (dim == 3) GmfGetLin(file, GmfVertices, &tmp[0], &tmp[1], &tmp[2], &ref);
    std::cout << "vert " << i << " ref " << ref << '\n';
    for (int j = 0; j < dim; ++j) coords[i * dim + j] = tmp[j];
  }
  GmfCloseMesh(file);
  (void)m;
}

static void read_meshb(osh::Mesh* m, const char* filepath) {
  int ver, dim;
  auto file = GmfOpenMesh(filepath, GmfRead, &ver, &dim);
  if (!file) {
    Omega_h_fail("could not open Meshb file %s\n", filepath);
  }
  OMEGA_H_CHECK(dim == 2 || dim == 3);
  switch (ver) {
    case 1: read_meshb_version<1>(m, file, ver, dim); return;
    case 2: read_meshb_version<2>(m, file, ver, dim); return;
    case 3: read_meshb_version<3>(m, file, ver, dim); return;
    case 4: read_meshb_version<4>(m, file, ver, dim); return;
  }
  Omega_h_fail("unknown libMeshb version %d\n", ver);
}

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  if (argc != 3) {
    std::cout << "usage: " << argv[0] << " input.meshb output.osh\n";
    return -1;
  }
  osh::Mesh m;
  read_meshb(&m, argv[1]);
  return 0;
}
