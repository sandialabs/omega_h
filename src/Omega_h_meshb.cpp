#include <libmeshb7.h>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
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
  HostWrite<LO> eqs2verts[4];
  HostWrite<ClassId> eqs2class_id[4];
  LO neqs[4];
  for (Int ent_dim = 1; ent_dim <= dim; ++ent_dim) {
    auto ent_kwd = simplex_kwds[ent_dim];
    neqs[ent_dim] = LO(GmfStatKwd(file, ent_kwd));
    if (ent_dim < dim && neqs[ent_dim] < 1) continue;
    eqs2verts[ent_dim] = HostWrite<LO>(LO(neqs[ent_dim]) * (ent_dim + 1));
    eqs2class_id[ent_dim] = HostWrite<ClassId>(LO(neqs[ent_dim]));
    safe_goto(file, ent_kwd);
    bool is_old_convention = true;
    for (LO i = 0; i < neqs[ent_dim]; ++i) {
      GmfIndex class_id;
      Few<GmfIndex, 4> tmp;
      if (ent_dim == 1) {
        GmfGetLin(file, ent_kwd, &tmp[0], &tmp[1], &class_id);
      } else if (ent_dim == 2) {
        GmfGetLin(file, ent_kwd, &tmp[0], &tmp[1], &tmp[2], &class_id);
      } else {
        GmfGetLin(file, ent_kwd, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &class_id);
      }
      for (Int j = 0; j < (ent_dim + 1); ++j) {
        eqs2verts[ent_dim][i * (ent_dim + 1) + j] = LO(tmp[j]) - 1;
      }
      eqs2class_id[ent_dim][i] = ClassId(class_id);
      if (LO(class_id) != i + 1) is_old_convention = false;
    }
    if (is_old_convention) {
      eqs2class_id[ent_dim] = HostWrite<ClassId>();
    }
  }
  GmfCloseMesh(file);
  build_from_elems2verts(
      mesh, OMEGA_H_SIMPLEX, dim, LOs(eqs2verts[dim].write()), LO(nverts));
  mesh->add_tag(VERT, "coordinates", dim, Reals(coords.write()));
  if (eqs2class_id[dim].exists()) {
    mesh->add_tag(dim, "class_id", 1, Read<ClassId>(ClassId(neqs[dim]), 1));
  } else {
    mesh->add_tag(dim, "class_id", 1, Read<ClassId>(eqs2class_id[dim].write()));
  }
  for (Int ent_dim = 1; ent_dim < dim; ++ent_dim) {
    if (eqs2class_id[ent_dim].exists()) {
      classify_equal_order(mesh, ent_dim, eqs2verts[ent_dim].write(),
          eqs2class_id[ent_dim].write());
    }
  }
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
    vert_refs = mesh->get_array<ClassId>(VERT, "class_id");
  } else {
    vert_refs = LOs(mesh->nverts(), 1);
  }
  auto h_vert_refs = HostRead<LO>(vert_refs);
  for (LO i = 0; i < nverts; ++i) {
    Few<GmfReal, 3> tmp = {0.0, 0.0, 0.0};
    for (int j = 0; j < dim; ++j) tmp[j] = GmfReal(coords[i * dim + j]);
    auto ref = GmfIndex(h_vert_refs[i]);
    if (dim == 2)
      GmfSetLin(file, GmfVertices, tmp[0], tmp[1], ref);
    else
      GmfSetLin(file, GmfVertices, tmp[0], tmp[1], tmp[2], ref);
  }
  for (Int ent_dim = 1; ent_dim <= dim; ++ent_dim) {
    auto ents2class_dim = mesh->get_array<I8>(ent_dim, "class_dim");
    auto ents2class_id = mesh->get_array<ClassId>(ent_dim, "class_id");
    auto ents2verts = mesh->ask_verts_of(ent_dim);
    auto ents_are_eqs = each_eq_to(ents2class_dim, I8(ent_dim));
    auto eqs2ents = collect_marked(ents_are_eqs);
    auto neqs = eqs2ents.size();
    auto eqs2verts = HostRead<LO>(unmap(eqs2ents, ents2verts, ent_dim + 1));
    auto eqs2class_id = HostRead<LO>(unmap(eqs2ents, ents2class_id, 1));
    auto ent_kwd = simplex_kwds[ent_dim];
    GmfSetKwd(file, ent_kwd, GmfLine(neqs));
    for (LO i = 0; i < neqs; ++i) {
      GmfIndex ref = eqs2class_id[i];
      Few<GmfIndex, 4> tmp;
      for (Int j = 0; j < (ent_dim + 1); ++j) {
        tmp[j] = eqs2verts[i * (ent_dim + 1) + j] + 1;
      }
      if (ent_dim == 1) {
        GmfSetLin(file, ent_kwd, tmp[0], tmp[1], ref);
      } else if (ent_dim == 2) {
        GmfSetLin(file, ent_kwd, tmp[0], tmp[1], tmp[2], ref);
      } else {
        GmfSetLin(file, ent_kwd, tmp[0], tmp[1], tmp[2], tmp[3], ref);
      }
    }
  }
  GmfCloseMesh(file);
}

void read(Mesh* mesh, std::string const& filepath) {
  int version, dim;
  auto file = GmfOpenMesh(filepath.c_str(), GmfRead, &version, &dim);
  if (!file) {
    Omega_h_fail(
        "could not open Meshb file %s for reading\n", filepath.c_str());
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

void write(Mesh* mesh, std::string const& filepath, int version) {
  auto dim = int(mesh->dim());
  auto file = GmfOpenMesh(filepath.c_str(), GmfWrite, version, dim);
  if (!file) {
    Omega_h_fail(
        "could not open Meshb file %s for writing\n", filepath.c_str());
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

template <int version>
static void read_sol_version(Mesh* mesh, GmfFile file, Int dim,
    std::string const& filepath, std::string const& sol_name) {
  using GmfReal = typename VersionTypes<version>::RealIn;
  int type_table[1];
  int ntypes, sol_size;
  LO nverts =
      LO(GmfStatKwd(file, GmfSolAtVertices, &ntypes, &sol_size, type_table));
  safe_goto(file, GmfSolAtVertices);
  if (ntypes != 1) {
    Omega_h_fail("\"%s\" has %d fields, Omega_h supports only one\n",
        filepath.c_str(), ntypes);
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
    Omega_h_fail(
        "unexpected field type %d in \"%s\"\n", field_type, filepath.c_str());
  }
  HostWrite<Real> hw(ncomps * nverts);
  for (LO i = 0; i < nverts; ++i) {
    Few<GmfReal, 6> tmp;
    GmfGetLin(file, GmfSolAtVertices, tmp.data());
    for (Int j = 0; j < ncomps; ++j) hw[i * ncomps + j] = Real(tmp[j]);
  }
  GmfCloseMesh(file);
  auto dr = Reals(hw.write());
  if (field_type == 3) dr = symms_inria2osh(dim, dr);
  mesh->add_tag(VERT, sol_name, ncomps, dr);
}

void read_sol(
    Mesh* mesh, std::string const& filepath, std::string const& sol_name) {
  int version, dim;
  auto file = GmfOpenMesh(filepath.c_str(), GmfRead, &version, &dim);
  if (!file) {
    Omega_h_fail("could not open Meshb solution file %s for reading\n",
        filepath.c_str());
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

template <int version>
static void write_sol_version(
    Mesh* mesh, GmfFile file, std::string const& sol_name) {
  auto dim = mesh->dim();
  using GmfReal = typename VersionTypes<version>::RealIn;
  auto nverts = mesh->nverts();
  auto tag = mesh->get_tag<Real>(VERT, sol_name);
  auto ncomps = tag->ncomps();
  int field_type = -1;
  if (ncomps == 1) {
    field_type = 1;
  } else if (ncomps == dim) {
    field_type = 2;
  } else if (ncomps == symm_ncomps(dim)) {
    field_type = 3;
  } else {
    Omega_h_fail(
        "unexpected # of components %d in tag %s\n", ncomps, sol_name.c_str());
  }
  int type_table[1] = {field_type};
  constexpr int ntypes = 1;
  GmfSetKwd(file, GmfSolAtVertices, GmfLine(nverts), ntypes, type_table);
  auto dr = tag->array();
  if (field_type == 3) dr = symms_osh2inria(dim, dr);
  HostRead<Real> hr(dr);
  for (LO i = 0; i < nverts; ++i) {
    Few<GmfReal, 6> tmp;
    for (Int j = 0; j < ncomps; ++j) tmp[j] = GmfReal(hr[i * ncomps + j]);
    GmfSetLin(file, GmfSolAtVertices, tmp.data());
  }
  GmfCloseMesh(file);
}

void write_sol(Mesh* mesh, std::string const& filepath,
    std::string const& sol_name, int version) {
  auto dim = int(mesh->dim());
  auto file = GmfOpenMesh(filepath.c_str(), GmfWrite, version, dim);
  if (!file) {
    Omega_h_fail("could not open Meshb solution file %s for writing\n",
        filepath.c_str());
  }
  OMEGA_H_CHECK(dim == 2 || dim == 3);
  switch (version) {
    case 1:
      write_sol_version<1>(mesh, file, sol_name);
      return;
    case 2:
      write_sol_version<2>(mesh, file, sol_name);
      return;
    case 3:
      write_sol_version<3>(mesh, file, sol_name);
      return;
    case 4:
      write_sol_version<4>(mesh, file, sol_name);
      return;
  }
  Omega_h_fail("unknown libMeshb solution version %d when reading\n", version);
}

}  // namespace meshb

}  // namespace Omega_h
