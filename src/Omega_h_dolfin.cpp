#include <Omega_h_dolfin.hpp>

#include <vector>

#include <Omega_h_adj.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_scan.hpp>
#include <Omega_h_shape.hpp>

namespace Omega_h {

static void form_sharing(Mesh* mesh_osh, Int ent_dim,
    std::map<std::int32_t, std::set<unsigned int>>* shared_ents) {
  auto n = mesh_osh->nents(ent_dim);
  if (!mesh_osh->could_be_shared(ent_dim)) {
    return;
  }
  auto dist = mesh_osh->ask_dist(ent_dim).invert();
  auto d_owners2copies = dist.roots2items();
  auto d_copies2rank = dist.items2ranks();
  auto d_copies2indices = dist.items2dest_idxs();
  auto h_owners2copies = HostRead<LO>(d_owners2copies);
  auto h_copies2rank = HostRead<I32>(d_copies2rank);
  auto h_copies2indices = HostRead<LO>(d_copies2indices);
  std::vector<I32> full_src_ranks;
  std::vector<I32> full_dest_ranks;
  std::vector<LO> full_dest_indices;
  auto my_rank = mesh_osh->comm()->rank();
  for (LO i_osh = 0; i_osh < n; ++i_osh) {
    auto begin = h_owners2copies[i_osh];
    auto end = h_owners2copies[i_osh + 1];
    if (end - begin <= 1) continue;
    for (LO copy = begin; copy < end; ++copy) {
      auto dest_rank = h_copies2rank[copy];
      auto dest_index = h_copies2indices[copy];
      for (LO copy2 = begin; copy2 < end; ++copy2) {
        auto src_rank = h_copies2rank[copy2];
        full_src_ranks.push_back(src_rank);
        full_dest_ranks.push_back(dest_rank);
        full_dest_indices.push_back(dest_index);
      }
    }
  }
  auto h_full_src_ranks = HostWrite<I32>(LO(full_src_ranks.size()));
  auto h_full_dest_ranks = HostWrite<I32>(LO(full_src_ranks.size()));
  auto h_full_dest_indices = HostWrite<I32>(LO(full_dest_indices.size()));
  for (LO i = 0; i < h_full_src_ranks.size(); ++i) {
    h_full_src_ranks[i] = full_src_ranks[size_t(i)];
    h_full_dest_ranks[i] = full_dest_ranks[size_t(i)];
    h_full_dest_indices[i] = full_dest_indices[size_t(i)];
  }
  auto d_full_src_ranks = Read<I32>(h_full_src_ranks.write());
  auto d_full_dest_ranks = Read<I32>(h_full_dest_ranks.write());
  auto d_full_dest_indices = Read<I32>(h_full_dest_indices.write());
  auto dist2 = Dist();
  dist2.set_parent_comm(mesh_osh->comm());
  dist2.set_dest_ranks(d_full_dest_ranks);
  dist2.set_dest_idxs(d_full_dest_indices, n);
  auto d_exchd_full_src_ranks = dist2.exch(d_full_src_ranks, 1);
  auto d_shared2ranks = dist2.invert().roots2items();
  auto h_exchd_full_src_ranks = HostRead<I32>(d_exchd_full_src_ranks);
  auto h_shared2ranks = HostRead<LO>(d_shared2ranks);
  auto ents_are_shared_w = HostWrite<Byte>(n);
  for (LO i_osh = 0; i_osh < n; ++i_osh) {
    auto begin = h_shared2ranks[i_osh];
    auto end = h_shared2ranks[i_osh + 1];
    ents_are_shared_w[i_osh] = ((end - begin) > 0);
  }
  for (LO i_osh = 0; i_osh < n; ++i_osh) {
    auto begin = h_shared2ranks[i_osh];
    auto end = h_shared2ranks[i_osh + 1];
    for (auto j = begin; j < end; ++j) {
      auto rank = h_exchd_full_src_ranks[j];
      if (rank != my_rank) {
        (*shared_ents)[i_osh].insert(unsigned(rank));
      }
    }
  }
}

void to_dolfin(dolfin::Mesh& mesh_dolfin, Mesh* mesh_osh) {
  dolfin::MeshEditor editor;
  static char const* const cell_type_names[4] = {
      "point", "interval", "triangle", "tetrahedron"};
  auto dim = mesh_osh->dim();
  OMEGA_H_CHECK(mesh_osh->parting() == OMEGA_H_ELEM_BASED);
  std::map<std::int32_t, std::set<unsigned int>> shared_verts;
  form_sharing(mesh_osh, VERT, &shared_verts);
  editor.open(mesh_dolfin, cell_type_names[dim], dim, dim);
  auto nverts = mesh_osh->nverts();
  auto nverts_global = mesh_osh->nglobal_ents(VERT);
  editor.init_vertices_global(nverts, nverts_global);
  auto d_vert_globals = mesh_osh->globals(VERT);
  auto h_vert_globals = HostRead<GO>(d_vert_globals);
  auto d_coords = mesh_osh->coords();
  auto h_coords = HostRead<Real>(d_coords);
  for (LO i_osh = 0; i_osh < nverts; ++i_osh) {
    editor.add_vertex_global(i_osh, h_vert_globals[i_osh],
        dolfin::Point(dim, &h_coords[i_osh * dim]));
  }
  auto ncells = mesh_osh->nelems();
  auto ncells_global = mesh_osh->nglobal_ents(dim);
  editor.init_cells_global(ncells, ncells_global);
  auto d_cell_globals = mesh_osh->globals(dim);
  auto h_cell_globals = HostRead<GO>(d_cell_globals);
  auto d_conn = mesh_osh->ask_elem_verts();
  auto h_conn = HostRead<LO>(d_conn);
  auto nverts_per_cell = dim + 1;
  /* this type needs to have .size(), .begin(), and .end().
     could have templated on dimension, but there is no threading
     here so its okay to just allocate it once. */
  std::vector<LO> cell_verts(nverts_per_cell);
  for (LO i = 0; i < ncells; ++i) {
    for (LO j = 0; j < (nverts_per_cell); ++j) {
      auto vert_osh = h_conn[i * (nverts_per_cell) + j];
      cell_verts[j] = vert_osh;
    }
    editor.add_cell(i, h_cell_globals[i], cell_verts);
  }
  /* the reordering in question is local to each cell,
     and among other things it requires the vertices of
     a cell to be ordered by increasing global number.
     this seems to be called "UFC order" */
  bool should_reorder = true;
  editor.close(should_reorder);

  // Set the ghost cell offset
  mesh_dolfin.topology().init_ghost(dim, size_t(ncells));

  // Set the ghost vertex offset
  mesh_dolfin.topology().init_ghost(0, size_t(nverts));

  // Assign map of shared cells and vertices
  mesh_dolfin.topology().shared_entities(0) = shared_verts;

  // Initialise number of globally connected cells to each facet. This
  // is necessary to distinguish between facets on an exterior
  // boundary and facets on a partition boundary (see
  // https://bugs.launchpad.net/dolfin/+bug/733834).
  dolfin::DistributedMeshTools::init_facet_cell_connections(mesh_dolfin);
}

/* due to their "UFC order" which requires vertices of a cell
   to be in increasing order of global ID,
   DOLFIN elements may look inverted from the Omega_h perspective.
   we can only assume that the DOLFIN connectivity is reasonable
   and adjust our own cell ordering to ensure positive volumes */
template <Int dim>
static void fix_inverted_elements_dim(Write<LO> elem_verts, Reals coords) {
  constexpr auto nverts_per_cell = dim + 1;
  auto nelems = divide_no_remainder(elem_verts.size(), nverts_per_cell);
  auto f = OMEGA_H_LAMBDA(LO e) {
    auto eev2v = gather_verts<nverts_per_cell>(elem_verts, e);
    auto eev2x = gather_vectors<nverts_per_cell, dim>(coords, eev2v);
    auto b = simplex_basis<dim, dim>(eev2x);
    auto s = simplex_size_from_basis(b);
    if (s < 0.0) {
      swap2(elem_verts[e * (nverts_per_cell) + 0],
          elem_verts[e * (nverts_per_cell) + 1]);
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
  auto nverts_per_cell = dim + 1;
  if (nelems)
    OMEGA_H_CHECK(Int(elem_verts_dolfin.size(0)) == (nverts_per_cell));
  auto h_elem_verts = HostWrite<LO>(nelems * (nverts_per_cell));
  for (LO i = 0; i < nelems; ++i) {
    auto ptr_dolfin = elem_verts_dolfin(i);
    for (Int j = 0; j < (nverts_per_cell); ++j) {
      h_elem_verts[i * (nverts_per_cell) + j] = ptr_dolfin[j];
    }
  }
  auto d_elem_verts = h_elem_verts.write();
  fix_inverted_elements(dim, d_elem_verts, d_coords);
  build_from_elems2verts(mesh_osh, OMEGA_H_SIMPLEX, dim, d_elem_verts, nverts);
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
    auto retrieve_indices = std::vector<std::size_t>(nents);
    for (LO i = 0; i < nents; ++i) retrieve_indices[i] = i;
    auto dof_indices =
        dofmap->entity_dofs(*mesh_dolfin, ent_dim, retrieve_indices);
    auto h_data = HostWrite<Real>(nents * ndofs_per_ent);
    vector->get_local(h_data.data(), nents * ndofs_per_ent, dof_indices.data());
    auto d_data = Reals(h_data.write());
    mesh_osh->remove_tag(ent_dim, name);
    mesh_osh->add_tag(ent_dim, name, ndofs_per_ent, d_data);
  }
}

}  // namespace Omega_h
