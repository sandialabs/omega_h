#include <libmeshb7.h>
#include "Omega_h.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_map.hpp"

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
static void read_meshb_version(Mesh* mesh, GmfFile file, int dim) {
  using GmfIndex = typename VersionTypes<version>::Index;
  using GmfReal = typename VersionTypes<version>::RealIn;
  LO nverts = LO(GmfStatKwd(file, GmfVertices));
  safe_goto(file, GmfVertices);
  auto coords = HostWrite<Real>(LO(nverts) * dim);
  for (LO i = 0; i < nverts; ++i) {
    GmfIndex ref;
    Few<GmfReal, 3> tmp;
    if (dim == 2)
      GmfGetLin(file, GmfVertices, &tmp[0], &tmp[1], &ref);
    else
      GmfGetLin(file, GmfVertices, &tmp[0], &tmp[1], &tmp[2], &ref);
    for (int j = 0; j < dim; ++j) coords[i * dim + j] = Real(tmp[j]);
  }
  auto side_kwd = simplex_kwds[dim - 1];
  LO nsides = LO(GmfStatKwd(file, side_kwd));
  auto sides2verts = HostWrite<LO>(LO(nsides) * dim);
  auto sides2class_id = HostWrite<I32>(LO(nsides));
  safe_goto(file, side_kwd);
  for (LO i = 0; i < nsides; ++i) {
    GmfIndex class_id;
    Few<GmfIndex, 3> tmp;
    if (dim == 2)
      GmfGetLin(file, side_kwd, &tmp[0], &tmp[1], &class_id);
    else
      GmfGetLin(file, side_kwd, &tmp[0], &tmp[1], &tmp[2], &class_id);
    for (int j = 0; j < dim; ++j) {
      sides2verts[i * dim + j] = LO(tmp[j] - 1);
    }
    sides2class_id[i] = LO(class_id);
  }
  auto elem_kwd = simplex_kwds[dim];
  LO nelems = LO(GmfStatKwd(file, elem_kwd));
  safe_goto(file, elem_kwd);
  auto elems2verts = HostWrite<LO>(LO(nelems) * (dim + 1));
  auto elem_class_ids = HostWrite<LO>(LO(nelems));
  auto is_old_convention = true;
  for (LO i = 0; i < nelems; ++i) {
    GmfIndex ref;
    Few<GmfIndex, 4> tmp;
    if (dim == 2)
      GmfGetLin(file, elem_kwd, &tmp[0], &tmp[1], &tmp[2], &ref);
    else
      GmfGetLin(file, elem_kwd, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &ref);
    if (ref != GmfIndex(i + 1)) is_old_convention = false;
    elem_class_ids[i] = LO(ref);
    for (int j = 0; j <= dim; ++j) {
      elems2verts[i * (dim + 1) + j] = LO(tmp[j] - 1);
    }
  }
  GmfCloseMesh(file);
  build_from_elems2verts(mesh, dim, LOs(elems2verts.write()), LO(nverts));
  mesh->add_tag(VERT, "coordinates", dim, Reals(coords.write()));
  if (is_old_convention) {
    mesh->add_tag(dim, "class_id", 1, LOs(LO(nelems), 1));
  } else {
    mesh->add_tag(dim, "class_id", 1, LOs(elem_class_ids.write()));
  }
  classify_equal_order(
      mesh, dim - 1, sides2verts.write(), sides2class_id.write());
  finalize_classification(mesh);
}

template <int version>
static void write_meshb_version(Mesh* mesh, GmfFile file, int dim) {
  using GmfIndex = typename VersionTypes<version>::Index;
  using GmfReal = typename VersionTypes<version>::RealOut;
  auto nverts = mesh->nverts();
  GmfSetKwd(file, GmfVertices, GmfLine(nverts));
  auto coords = HostRead<Real>(mesh->coords());
  LOs vert_refs;
  if (mesh->has_tag(VERT, "class_id")) {
    vert_refs = mesh->get_array<LO>(VERT, "class_id");
  } else {
    vert_refs = LOs(mesh->nverts(), 1);
  }
  auto h_vert_refs = HostRead<LO>(vert_refs);
  for (LO i = 0; i < nverts; ++i) {
    Few<GmfReal, 3> tmp;
    for (int j = 0; j < dim; ++j) tmp[j] = GmfReal(coords[i * dim + j]);
    auto ref = GmfIndex(h_vert_refs[i]);
    if (dim == 2)
      GmfSetLin(file, GmfVertices, tmp[0], tmp[1], ref);
    else
      GmfSetLin(file, GmfVertices, tmp[0], tmp[1], tmp[2], ref);
  }
  auto side_kwd = simplex_kwds[dim - 1];
  auto sds2class_dim = mesh->get_array<I8>(dim - 1, "class_dim");
  auto sds_are_sides = each_eq_to(sds2class_dim, I8(dim - 1));
  auto sides2sds = collect_marked(sds_are_sides);
  auto nsides = sides2sds.size();
  auto sds2verts = mesh->ask_verts_of(dim - 1);
  auto sds2class_id = mesh->get_array<LO>(dim - 1, "class_id");
  auto sides2verts = HostRead<LO>(unmap(sides2sds, sds2verts, dim));
  auto sides2class_id = HostRead<LO>(unmap(sides2sds, sds2class_id, 1));
  GmfSetKwd(file, side_kwd, GmfLine(nsides));
  for (LO i = 0; i < nsides; ++i) {
    GmfIndex ref = sides2class_id[i];
    Few<GmfIndex, 3> tmp;
    for (int j = 0; j < dim; ++j) tmp[j] = sides2verts[i * dim + j] + 1;
    if (dim == 2)
      GmfSetLin(file, side_kwd, tmp[0], tmp[1], ref);
    else
      GmfSetLin(file, side_kwd, tmp[0], tmp[1], tmp[2], ref);
  }
  auto elem_kwd = simplex_kwds[dim];
  auto nelems = mesh->nelems();
  auto elems2verts = HostRead<LO>(mesh->ask_elem_verts());
  auto elem_class_ids = HostRead<LO>(mesh->get_array<LO>(dim, "class_id"));
  GmfSetKwd(file, elem_kwd, GmfLine(nelems));
  for (LO i = 0; i < nelems; ++i) {
    auto ref = GmfIndex(elem_class_ids[i]);
    Few<GmfIndex, 4> tmp;
    for (int j = 0; j < dim + 1; ++j) {
      tmp[j] = GmfIndex(elems2verts[i * (dim + 1) + j] + 1);
    }
    if (dim == 2)
      GmfSetLin(file, elem_kwd, tmp[0], tmp[1], tmp[2], ref);
    else
      GmfSetLin(file, elem_kwd, tmp[0], tmp[1], tmp[2], tmp[3], ref);
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

OMEGA_H_INLINE static Matrix<1, 1> vector2symm_inria(Vector<1> v) {
  return matrix_1x1(v[0]);
}

OMEGA_H_INLINE static Matrix<2, 2> vector2symm_inria(Vector<3> v) {
  Matrix<2, 2> symm;
  symm[0][0] = v[0];
  symm[0][1] = v[1];
  symm[1][1] = v[2];
  symm[1][0] = symm[0][1];
  return symm;
}

OMEGA_H_INLINE static Matrix<3, 3> vector2symm_inria(Vector<6> v) {
  Matrix<3, 3> symm;
  symm[0][0] = v[0];
  symm[0][1] = v[1];
  symm[1][1] = v[2];
  symm[0][2] = v[3];
  symm[1][2] = v[4];
  symm[2][2] = v[5];
  symm[1][0] = symm[0][1];
  symm[2][0] = symm[0][2];
  symm[2][1] = symm[1][2];
  return symm;
}

template <Int dim>
static Reals symms_inria2osh_dim(Reals symms) {
  auto n = divide_no_remainder(symms.size(), symm_ncomps(dim));
  Write<Real> out(symms.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto iv = get_vector<symm_ncomps(dim)>(symms, i);
    auto is = vector2symm_inria(iv);
    auto ov = symm2vector(is);
    set_vector(out, i, ov);
  };
  parallel_for(n, f);
  return out;
}

static Reals symms_inria2osh(Int dim, Reals symms) {
  if (dim == 3) return symms_inria2osh_dim<3>(symms);
  if (dim == 2) return symms_inria2osh_dim<2>(symms);
  if (dim == 1) return symms_inria2osh_dim<1>(symms);
  OMEGA_H_NORETURN(Reals());
}

#if 0
/* Symmetric metric tensor storage convention used by
   INRIA:  https://hal.inria.fr/inria-00363007/document */
OMEGA_H_INLINE static Vector<1> symm2vector_inria(Matrix<1, 1> symm) {
  return vector_1(symm[0][0]);
}

OMEGA_H_INLINE static Vector<3> symm2vector_inria(Matrix<2, 2> symm) {
  Vector<3> v;
  v[0] = symm[0][0];
  v[1] = symm[0][1];
  v[2] = symm[1][1];
  return v;
}

OMEGA_H_INLINE static Vector<6> symm2vector_inria(Matrix<3, 3> symm) {
  Vector<6> v;
  v[0] = symm[0][0];
  v[1] = symm[0][1];
  v[2] = symm[1][1];
  v[3] = symm[0][2];
  v[4] = symm[1][2];
  v[5] = symm[2][2];
  return v;
}

template <Int dim>
static Reals symms_osh2inria_dim(Reals symms) {
  auto n = divide_no_remainder(symms.size(), symm_ncomps(dim));
  Write<Real> out(symms.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto iv = get_vector<symm_ncomps(dim)>(symms, i);
    auto is = vector2symm(iv);
    auto ov = symm2vector_inria(is);
    set_vector(out, i, ov);
  };
  parallel_for(n, f);
  return out;
}

static Reals symms_osh2inria(Int dim, Reals symms) {
  if (dim == 3) return symms_osh2inria_dim<3>(symms);
  if (dim == 2) return symms_osh2inria_dim<2>(symms);
  if (dim == 1) return symms_osh2inria_dim<1>(symms);
  OMEGA_H_NORETURN(Reals());
}
#endif

template <int version>
static void read_sol_version(Mesh* mesh, GmfFile file, Int dim,
    const char* filepath, const char* sol_name) {
  using GmfReal = typename VersionTypes<version>::RealIn;
  int type_table[1];
  int ntypes, sol_size;
  LO nverts =
      LO(GmfStatKwd(file, GmfSolAtVertices, &ntypes, &sol_size, type_table));
  safe_goto(file, GmfSolAtVertices);
  if (ntypes != 1) {
    Omega_h_fail(
        "\"%s\" has %d fields, Omega_h supports only one\n", filepath, ntypes);
    /* also, there was definitely a buffer overflow writing to type_table */
  }
  auto field_type = type_table[0];
  Int ncomps = -1;
  if (field_type == 1)
    ncomps = 1;
  else if (field_type == 2)
    ncomps = dim;
  else if (field_type == 3)
    ncomps = symm_ncomps(dim);
  else {
    Omega_h_fail("unexpected field type %d in \"%s\"\n", field_type, filepath);
  }
  HostWrite<Real> hw(ncomps * nverts);
  for (LO i = 0; i < nverts; ++i) {
    Few<GmfReal, 6> tmp;
    if (ncomps == 1)
      GmfGetLin(file, GmfSolAtVertices, &tmp[0]);
    else if (ncomps == 2)
      GmfGetLin(file, GmfSolAtVertices, &tmp[0], &tmp[1]);
    else if (ncomps == 3)
      GmfGetLin(file, GmfSolAtVertices, &tmp[0], &tmp[1], &tmp[2]);
    else if (ncomps == 6) {
      GmfGetLin(file, GmfSolAtVertices, &tmp[0], &tmp[1], &tmp[2], &tmp[3],
          &tmp[4], &tmp[5]);
    }
    for (Int j = 0; j < ncomps; ++j) hw[i * ncomps + j] = Real(tmp[j]);
  }
  GmfCloseMesh(file);
  auto dr = Reals(hw.write());
  if (field_type == 3) dr = symms_inria2osh(dim, dr);
  mesh->add_tag(VERT, sol_name, ncomps, dr);
}

void read_sol(Mesh* mesh, const char* filepath, const char* sol_name) {
  int version, dim;
  auto file = GmfOpenMesh(filepath, GmfRead, &version, &dim);
  if (!file) {
    Omega_h_fail(
        "could not open Meshb solution file %s for reading\n", filepath);
  }
  OMEGA_H_CHECK(dim == 2 || dim == 3);
  switch (version) {
    case 1:
      read_sol_version<1>(mesh, file, dim, filepath, sol_name);
      return;
    case 2:
      read_sol_version<2>(mesh, file, dim, filepath, sol_name);
      return;
    case 3:
      read_sol_version<3>(mesh, file, dim, filepath, sol_name);
      return;
    case 4:
      read_sol_version<4>(mesh, file, dim, filepath, sol_name);
      return;
  }
  Omega_h_fail("unknown libMeshb solution version %d when reading\n", version);
}

}  // namespace meshb

}  // namespace Omega_h
