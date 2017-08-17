#include <Omega_h_dolfin.hpp>

#include <vector>

#include <Omega_h_adj.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_loop.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>

namespace Omega_h {

void to_dolfin(dolfin::Mesh& mesh_dolfin, Mesh* mesh_osh) {
  dolfin::MeshEditor editor;
  static char const* const cell_type_names[4] = {
      "point", "interval", "triangle", "tetrahedron"};
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
    editor.add_vertex_global(
        i, h_vert_globals[i], dolfin::Point(dim, &h_coords[i * dim]));
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
  std::vector<LO> cell_verts(dim + 1);
  for (LO i = 0; i < ncells; ++i) {
    for (LO j = 0; j < (dim + 1); ++j) {
      cell_verts[j] = h_conn[i * (dim + 1) + j];
    }
    editor.add_cell(i, h_cell_globals[i], cell_verts);
  }
  bool should_reorder = false;
  editor.close(should_reorder);
}

/* DOLFIN intermixes inverted elements with non-inverted ones!
   best we can do is to reverse the ordering of
   the inverted ones */
template <Int dim>
static void fix_inverted_elements_dim(Write<LO> elem_verts, Reals coords) {
  auto nelems = divide_no_remainder(elem_verts.size(), dim + 1);
  auto f = OMEGA_H_LAMBDA(LO e) {
    auto eev2v = gather_verts<dim + 1>(elem_verts, e);
    auto eev2x = gather_vectors<dim + 1, dim>(coords, eev2v);
    auto b = simplex_basis<dim, dim>(eev2x);
    auto s = element_size(b);
    if (s < 0.0) {
      swap2(elem_verts[e * (dim + 1) + 0], elem_verts[e * (dim + 1) + 1]);
    }
  };
  parallel_for(nelems, f, "fix_inverted_elements");
}

static void fix_inverted_elements(Int dim, Write<LO> elem_verts, Reals coords) {
  if (dim == 1) fix_inverted_elements_dim<1>(elem_verts, coords);
  if (dim == 2) fix_inverted_elements_dim<2>(elem_verts, coords);
  if (dim == 3) fix_inverted_elements_dim<3>(elem_verts, coords);
}

void from_dolfin(Mesh* mesh_osh, dolfin::Mesh const& mesh_dolfin) {
  auto& topology = mesh_dolfin.topology();
  auto& geometry = mesh_dolfin.geometry();
  OMEGA_H_CHECK(geometry.degree() == 1);
  OMEGA_H_CHECK(topology.dim() == geometry.dim());
  auto dim = Int(topology.dim());
  auto nverts = LO(topology.size(VERT));
  auto nelems = LO(topology.size(dim));
  auto& coords_dolfin = geometry.x();
  auto h_coords = HostWrite<Real>(nverts * dim);
  for (LO i = 0; i < nverts * dim; ++i) {
    h_coords[i] = coords_dolfin[i];
  }
  auto d_coords = Reals(h_coords.write());
  auto& vert_globals_dolfin = topology.global_indices(VERT);
  auto h_vert_globals = HostWrite<GO>(nverts);
  for (LO i = 0; i < nverts; ++i) {
    h_vert_globals[i] = vert_globals_dolfin[i];
  }
  auto d_vert_globals = GOs(h_vert_globals.write());
  auto& elem_verts_dolfin = topology(dim, VERT);
  if (nelems) OMEGA_H_CHECK(Int(elem_verts_dolfin.size(0)) == (dim + 1));
  auto h_elem_verts = HostWrite<LO>(nelems * (dim + 1));
  for (LO i = 0; i < nelems; ++i) {
    auto ptr_dolfin = elem_verts_dolfin(i);
    for (Int j = 0; j < (dim + 1); ++j) {
      h_elem_verts[i * (dim + 1) + j] = ptr_dolfin[j];
    }
  }
  auto d_elem_verts = h_elem_verts.write();
  fix_inverted_elements(dim, d_elem_verts, d_coords);
  build_from_elems2verts(mesh_osh, dim, d_elem_verts, nverts);
  mesh_osh->remove_tag(VERT, "global");
  mesh_osh->add_tag(VERT, "global", 1, d_vert_globals);
  mesh_osh->add_tag(VERT, "coordinates", dim, d_coords);
}

void from_dolfin(
    Mesh* mesh_osh, dolfin::Function const& function, std::string const& name) {
  auto function_space = function.function_space();
  auto vector = function.vector();
  auto dofmap = function_space->dofmap();
  auto mesh_dolfin = function_space->mesh();
  for (Int ent_dim = 0; ent_dim <= mesh_osh->dim(); ++ent_dim) {
    auto ndofs_per_ent = dofmap->num_entity_dofs(ent_dim);
    if (ndofs_per_ent == 0) continue;
    auto nents = mesh_osh->nents(ent_dim);
    auto entity_indices = std::vector<std::size_t>(nents);
    for (LO i = 0; i < nents; ++i) entity_indices[i] = i;
    auto dof_indices =
        dofmap->entity_dofs(*mesh_dolfin, ent_dim, entity_indices);
    auto h_data = HostWrite<Real>(nents * ndofs_per_ent);
    vector->get_local(h_data.data(), nents * ndofs_per_ent, dof_indices.data());
    auto d_data = Reals(h_data.write());
    mesh_osh->remove_tag(ent_dim, name);
    mesh_osh->add_tag(ent_dim, name, ndofs_per_ent, d_data);
  }
}

}  // namespace Omega_h
