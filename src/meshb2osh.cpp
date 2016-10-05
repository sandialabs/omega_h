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

static void safe_goto(long long file, GmfKwdCod key) {
  auto ret = GmfGotoKwd(file, key);
  OMEGA_H_CHECK(ret);
}

template <int version>
static void read_meshb_version(osh::Library const& lib,
    osh::Mesh* mesh, long long file, int ver, int dim) {
  using Index = typename VersionTypes<version>::Index;
  using Real = typename VersionTypes<version>::Real;
  using Line = typename VersionTypes<version>::Line;
  OMEGA_H_CHECK(ver == version);
  Index nverts = GmfStatKwd(file, GmfVertices);
  static GmfKwdCod const simplex_kwds[4] = {
    GmfVertices,
    GmfEdges,
    GmfTriangles,
    GmfTetrahedra
  };
  safe_goto(file, GmfVertices);
  auto coords = osh::HostWrite<osh::Real>(nverts * dim);
  for (Line i = 0; i < nverts; ++i) {
    Index ref;
    osh::Few<Real, 3> tmp;
    if (dim == 2) GmfGetLin(file, GmfVertices, &tmp[0], &tmp[1], &ref);
    if (dim == 3) GmfGetLin(file, GmfVertices, &tmp[0], &tmp[1], &tmp[2], &ref);
    OMEGA_H_CHECK(ref == i + 1);
    for (int j = 0; j < dim; ++j) coords[i * dim + j] = tmp[j];
  }
  auto side_kwd = simplex_kwds[dim - 1];
  Index nsides = GmfStatKwd(file, side_kwd);
  auto sides2verts = osh::HostWrite<osh::LO>(nsides * dim);
  auto sides2class_id = osh::HostWrite<osh::I32>(nsides);
  if (nsides) {
    safe_goto(file, side_kwd);
    for (Line i = 0; i < nsides; ++i) {
      Index class_id;
      osh::Few<Index, 3> tmp;
      if (dim == 2) GmfGetLin(file, side_kwd, &tmp[0], &tmp[1], &class_id);
      if (dim == 3) GmfGetLin(file, side_kwd, &tmp[0], &tmp[1], &tmp[2], &class_id);
      for (int j = 0; j < dim; ++j) sides2verts[i * dim + j] = tmp[j] - 1;
      sides2class_id[i] = class_id;
    }
  }
  auto elem_kwd = simplex_kwds[dim];
  Index nelems = GmfStatKwd(file, elem_kwd);
  safe_goto(file, elem_kwd);
  auto elems2verts = osh::HostWrite<osh::LO>(nelems * (dim + 1));
  for (Line i = 0; i < nelems; ++i) {
    Index ref;
    osh::Few<Index, 4> tmp;
    if (dim == 2) GmfGetLin(file, elem_kwd, &tmp[0], &tmp[1], &tmp[2], &ref);
    if (dim == 3) GmfGetLin(file, elem_kwd, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &ref);
    OMEGA_H_CHECK(ref == i + 1);
    for (int j = 0; j <= dim; ++j) elems2verts[i * (dim + 1) + j] = tmp[j] - 1;
  }
  osh::build_from_elems2verts(mesh, lib, dim,
      osh::LOs(elems2verts.write()), nverts);
  mesh->add_tag(osh::VERT, "coordinates", dim, OMEGA_H_LINEAR_INTERP, OMEGA_H_DO_OUTPUT,
      osh::Reals(coords.write()));
  GmfCloseMesh(file);
  (void)mesh;
}

static void read_meshb(osh::Library const& lib, osh::Mesh* mesh, const char* filepath) {
  int ver, dim;
  auto file = GmfOpenMesh(filepath, GmfRead, &ver, &dim);
  if (!file) {
    Omega_h_fail("could not open Meshb file %s\n", filepath);
  }
  OMEGA_H_CHECK(dim == 2 || dim == 3);
  switch (ver) {
    case 1: read_meshb_version<1>(lib, mesh, file, ver, dim); return;
    case 2: read_meshb_version<2>(lib, mesh, file, ver, dim); return;
    case 3: read_meshb_version<3>(lib, mesh, file, ver, dim); return;
    case 4: read_meshb_version<4>(lib, mesh, file, ver, dim); return;
  }
  Omega_h_fail("unknown libMeshb version %d\n", ver);
}

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  if (argc != 3) {
    std::cout << "usage: " << argv[0] << " input.meshb output.osh\n";
    return -1;
  }
  osh::Mesh mesh;
  read_meshb(lib, &mesh, argv[1]);
  osh::binary::write(argv[2], &mesh);
  return 0;
}
