#include <Omega_h_dolfin.hpp>

#include <vector>

#include <Omega_h_mesh.hpp>

#ifdef OMEGA_H_USE_MPI
#define HAS_MPI // omg DOLFIN this is so bad
#endif

#include <dolfin/mesh/MeshEditor.h>

namespace Omega_h {

void to_dolfin(dolfin::Mesh& mesh_dolfin, Mesh* mesh_osh) {
  dolfin::MeshEditor editor;
  static char const* const cell_type_names[4] = {
    "point",
    "interval",
    "triangle",
    "tetrahedron"
  };
  OMEGA_H_CHECK(mesh_osh->parting() == OMEGA_H_ELEM_BASED);
  auto dim = mesh_osh->dim();
  editor.open(mesh_dolfin, cell_type_names[dim], dim, dim);
  auto nverts = mesh_osh->nverts();
  auto nverts_global = mesh_osh->nglobal_ents(VERT);
  editor.init_vertices_global(nverts, nverts_global);
  auto d_vert_globals = mesh_osh->globals(VERT);
  auto h_vert_globals = HostRead<GO>(d_vert_globals);
  auto d_coords = mesh_osh->coords();
  auto h_coords = HostRead<Real>(d_coords);
  for (LO i = 0; i < nverts; ++i) {
    editor.add_vertex_global(i, h_vert_globals[i],
        dolfin::Point(dim, &h_coords[i * dim]));
  }
  auto ncells = mesh_osh->nelems();
  auto ncells_global = mesh_osh->nglobal_ents(dim);
  editor.init_cells_global(ncells, ncells_global);
  auto d_cell_globals = mesh_osh->globals(dim);
  auto h_cell_globals = HostRead<GO>(d_cell_globals);
  auto d_conn = mesh_osh->ask_elem_verts();
  auto h_conn = HostRead<LO>(d_conn);
  /* this type needs to have .size(), .begin(), and .end().
     could have templated on dimension, but there is no threading
     here so its okay to just allocate it once. */
  std::vector<LO> cell_verts(dim);
  for (LO i = 0; i < ncells; ++i) {
    for (LO j = 0; j < (dim + 1); ++j) {
      cell_verts[j] = h_conn[i * (dim + 1) + j];
    }
    editor.add_cell(i, h_cell_globals[i], cell_verts);
  }
  bool should_reorder = false;
  editor.close(should_reorder);
}

void from_dolfin(Mesh* mesh_osh, dolfin::Mesh const& mesh_dolfin) {
  (void)mesh_osh;
  (void)mesh_dolfin;
}

}
