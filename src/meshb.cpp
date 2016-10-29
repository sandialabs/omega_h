#include <Omega_h.hpp>
#include <libmeshb7.h>
#include "classify.hpp"
#include "map.hpp"

namespace Omega_h {

namespace meshb {

template <int version>
struct VersionTypes;

using GmfFile = std::int64_t;
/* contrary to the documentation, GmfSetKwd
 * and GmfStatKwd always operate on 64-bit values.
 */
using GmfLine = std::int64_t;

template <>
struct VersionTypes<1> {
  using Index = std::int32_t;
  using RealIn = float;
  // contrary to the documentation, GmfSetLin expects doubles
  using RealOut = double;
};

template <>
struct VersionTypes<2> {
  using Index = std::int32_t;
  using RealIn = double;
  using RealOut = double;
};

template <>
struct VersionTypes<3> {
  using Index = std::int32_t;
  using RealIn = double;
  using RealOut = double;
};

template <>
struct VersionTypes<4> {
  using Index = std::int64_t;
  using RealIn = double;
  using RealOut = double;
};

static void safe_goto(GmfFile file, GmfKwdCod key) {
  auto ret = GmfGotoKwd(file, key);
  OMEGA_H_CHECK(ret);
}

static GmfKwdCod const simplex_kwds[4] = {
    GmfVertices, GmfEdges, GmfTriangles, GmfTetrahedra};

template <int version>
static void read_meshb_version(
    Mesh* mesh, GmfFile file, int dim) {
  using GmfIndex = typename VersionTypes<version>::Index;
  using GmfReal = typename VersionTypes<version>::RealIn;
  auto nverts = GmfStatKwd(file, GmfVertices);
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
  auto nsides = GmfStatKwd(file, side_kwd);
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
  auto nelems = GmfStatKwd(file, elem_kwd);
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
  GmfCloseMesh(file);
  build_from_elems2verts(mesh, dim, LOs(elems2verts.write()), nverts);
  mesh->add_tag(VERT, "coordinates", dim, OMEGA_H_LINEAR_INTERP,
      OMEGA_H_DO_OUTPUT, Reals(coords.write()));
  classify_equal_order(mesh, dim - 1, sides2verts.write(), sides2class_id.write());
  finalize_classification(mesh);
}

template <int version>
static void write_meshb_version(
    Mesh* mesh, GmfFile file, int dim) {
  using GmfIndex = typename VersionTypes<version>::Index;
  using GmfReal = typename VersionTypes<version>::RealOut;
  GmfLine nverts = mesh->nverts();
  GmfSetKwd(file, GmfVertices, nverts);
  auto coords = HostRead<Real>(mesh->coords());
  for (GmfLine i = 0; i < nverts; ++i) {
    Few<GmfReal, 3> tmp;
    for (int j = 0; j < dim; ++j) tmp[j] = coords[i * dim + j];
    GmfIndex ref = i + 1;
    if (dim == 2) GmfSetLin(file, GmfVertices, tmp[0], tmp[1], ref);
    if (dim == 3) GmfSetLin(file, GmfVertices, tmp[0], tmp[1], tmp[2], ref);
  }
  auto side_kwd = simplex_kwds[dim - 1];
  auto sds2class_dim = mesh->get_array<I8>(dim - 1, "class_dim");
  auto sds_are_sides = each_eq_to(sds2class_dim, I8(dim - 1));
  auto sides2sds = collect_marked(sds_are_sides);
  GmfLine nsides = sides2sds.size();
  auto sds2verts = mesh->ask_verts_of(dim - 1);
  auto sds2class_id = mesh->get_array<LO>(dim - 1, "class_id");
  auto sides2verts = HostRead<LO>(unmap(sides2sds, sds2verts, dim));
  auto sides2class_id = HostRead<LO>(unmap(sides2sds, sds2class_id, 1));
  GmfSetKwd(file, side_kwd, nsides);
  for (GmfLine i = 0; i < nsides; ++i) {
    GmfIndex ref = sides2class_id[i];
    Few<GmfIndex, 3> tmp;
    for (int j = 0; j < dim; ++j) tmp[j] = sides2verts[i * dim + j];
    if (dim == 2) GmfSetLin(file, side_kwd, tmp[0], tmp[1], ref);
    if (dim == 3) GmfSetLin(file, side_kwd, tmp[0], tmp[1], tmp[2], ref);
  }
  auto elem_kwd = simplex_kwds[dim];
  GmfLine nelems = mesh->nelems();
  GmfSetKwd(file, elem_kwd, nelems);
  for (GmfLine i = 0; i < nelems; ++i) {
    GmfIndex ref = i + 1;
    Few<GmfIndex, 4> tmp;
    for (int j = 0; j < dim + 1; ++j) tmp[j] = sides2verts[i * (dim + 1) + j];
    if (dim == 2) GmfSetLin(file, elem_kwd, tmp[0], tmp[1], tmp[2], ref);
    if (dim == 3) GmfSetLin(file, elem_kwd, tmp[0], tmp[1], tmp[2], tmp[3], ref);
  }
  GmfCloseMesh(file);
}

void read(Mesh* mesh, const char* filepath) {
  int version, dim;
  auto file = GmfOpenMesh(filepath, GmfRead, &version, &dim);
  if (!file) {
    Omega_h_fail("could not open Meshb file %s for reading\n", filepath);
  }
  OMEGA_H_CHECK(dim == 2 || dim == 3);
  switch (version) {
    case 1:
      read_meshb_version<1>(mesh, file, dim);
      return;
    case 2:
      read_meshb_version<2>(mesh, file, dim);
      return;
    case 3:
      read_meshb_version<3>(mesh, file, dim);
      return;
    case 4:
      read_meshb_version<4>(mesh, file, dim);
      return;
  }
  Omega_h_fail("unknown libMeshb version %d when reading\n", version);
}

void write(Mesh* mesh, const char* filepath, int version) {
  auto dim = mesh->dim();
  auto file = GmfOpenMesh(filepath, GmfWrite, version, dim);
  if (!file) {
    Omega_h_fail("could not open Meshb file %s for writing\n", filepath);
  }
  switch (version) {
    case 1:
      write_meshb_version<1>(mesh, file, dim);
      return;
    case 2:
      write_meshb_version<2>(mesh, file, dim);
      return;
    case 3:
      write_meshb_version<3>(mesh, file, dim);
      return;
    case 4:
      write_meshb_version<4>(mesh, file, dim);
      return;
  }
  Omega_h_fail("unknown libMeshb version %d when writing\n", version);
}

}

}
