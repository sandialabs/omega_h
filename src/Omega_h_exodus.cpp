#include <exodusII.h>

#include <algorithm>
#include <set>

#include "Omega_h.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_internal.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

#define CALL(f) CHECK((f) >= 0)

namespace exodus {

static bool is_type_supported(int dim, std::string const& type) {
  if (dim == 2) {
    if (type == "tri3") return true;
    if (type == "TRI") return true;
  }
  if (dim == 3) {
    if (type == "tetra4") return true;
    if (type == "TETRA") return true;
  }
  return false;
}

// subtracts one and maps from Exodus
// side ordering to Omega_h
static INLINE int side_exo2osh(int dim, int side) {
  switch (dim) {
    case 2:
      switch (side) {
        case 1:
          return 0;
        case 2:
          return 1;
        case 3:
          return 2;
      }
    case 3:
      switch (side) {
        case 1:
          return 1;
        case 2:
          return 2;
        case 3:
          return 3;
        case 4:
          return 0;
      }
  }
  return -1;
}

// from Omega_h
// side ordering to Exodus and adds one
static INLINE int side_osh2exo(int dim, int side) {
  switch (dim) {
    case 2:
      switch (side) {
        case 0:
          return 1;
        case 1:
          return 2;
        case 2:
          return 3;
      }
    case 3:
      switch (side) {
        case 1:
          return 1;
        case 2:
          return 2;
        case 3:
          return 3;
        case 0:
          return 4;
      }
  }
  return -1;
}

void read(
    std::string const& path, Mesh* mesh, bool verbose, int classify_with) {
  auto comp_ws = int(sizeof(Real));
  int io_ws = 0;
  float version;
  auto mode = EX_READ | EX_MAPS_INT64_API;
  auto file = ex_open(path.c_str(), mode, &comp_ws, &io_ws, &version);
  if (file < 0) Omega_h_fail("can't open Exodus file %s\n", path.c_str());
  ex_init_params init_params;
  CALL(ex_get_init_ext(file, &init_params));
  if (verbose) {
    std::cout << "init params for " << path << ":\n";
    std::cout << " ExodusII " << version << '\n';
    std::cout << " Exodus ID " << file << '\n';
    std::cout << " comp_ws " << comp_ws << '\n';
    std::cout << " io_ws " << io_ws << '\n';
    std::cout << " Title " << init_params.title << '\n';
    std::cout << " num_dim " << init_params.num_dim << '\n';
    std::cout << " num_nodes " << init_params.num_nodes << '\n';
    std::cout << " num_elem " << init_params.num_elem << '\n';
    std::cout << " num_elem_blk " << init_params.num_elem_blk << '\n';
    std::cout << " num_node_sets " << init_params.num_node_sets << '\n';
    std::cout << " num_side_sets " << init_params.num_side_sets << '\n';
  }
  auto dim = int(init_params.num_dim);
  HostWrite<Real> h_coord_blk[3];
  for (Int i = 0; i < dim; ++i) {
    h_coord_blk[i] = HostWrite<Real>(init_params.num_nodes);
  }
  CALL(ex_get_coord(file, h_coord_blk[0].data(), h_coord_blk[1].data(),
      h_coord_blk[2].data()));
  HostWrite<Real> h_coords(init_params.num_nodes * dim);
  for (LO i = 0; i < init_params.num_nodes; ++i) {
    for (Int j = 0; j < dim; ++j) {
      h_coords[i * dim + j] = h_coord_blk[j][i];
    }
  }
  auto coords = Reals(h_coords.write());
  std::vector<int> block_ids(init_params.num_elem_blk);
  CALL(ex_get_ids(file, EX_ELEM_BLOCK, block_ids.data()));
  HostWrite<LO> h_conn(init_params.num_elem * (dim + 1));
  Write<LO> elem_class_ids_w(init_params.num_elem);
  LO start = 0;
  for (size_t i = 0; i < block_ids.size(); ++i) {
    char elem_type[MAX_STR_LENGTH + 1];
    int nentries;
    int nnodes_per_entry;
    int nedges_per_entry;
    int nfaces_per_entry;
    int nattr_per_entry;
    CALL(ex_get_block(file, EX_ELEM_BLOCK, block_ids[i], elem_type, &nentries,
        &nnodes_per_entry, &nedges_per_entry, &nfaces_per_entry,
        &nattr_per_entry));
    if (verbose) {
      std::cout << "block " << block_ids[i] << " has " << nentries
                << " elements of type " << elem_type << '\n';
    }
    if (!is_type_supported(dim, elem_type)) {
      Omega_h_fail("type %s is not supported for %dD !\n", elem_type, dim);
    }
    CHECK(nnodes_per_entry == dim + 1);
    if (nedges_per_entry < 0) nedges_per_entry = 0;
    if (nfaces_per_entry < 0) nfaces_per_entry = 0;
    std::vector<int> edge_conn(nentries * nedges_per_entry);
    std::vector<int> face_conn(nentries * nfaces_per_entry);
    CALL(ex_get_conn(file, EX_ELEM_BLOCK, block_ids[i], h_conn.data() + start,
        edge_conn.data(), face_conn.data()));
    auto region_id = block_ids[i];
    auto f0 = OMEGA_H_LAMBDA(LO entry) { elem_class_ids_w[start + entry] = region_id; };
    parallel_for(nentries, f0);
    start += nentries * nnodes_per_entry;
  }
  CHECK(start == init_params.num_elem * (dim + 1));
  auto conn = subtract_from_each(LOs(h_conn.write()), 1);
  build_from_elems_and_coords(mesh, dim, conn, coords);
  classify_elements(mesh);
  std::vector<int> side_set_ids(init_params.num_side_sets);
  CALL(ex_get_ids(file, EX_SIDE_SET, side_set_ids.data()));
  Write<LO> side_class_ids_w(mesh->nents(dim - 1), -1);
  auto sides_are_exposed = mark_exposed_sides(mesh);
  classify_sides_by_exposure(mesh, sides_are_exposed);
  Write<I8> side_class_dims_w =
      deep_copy(mesh->get_array<I8>(dim - 1, "class_dim"));
  auto exposed_sides2side = collect_marked(sides_are_exposed);
  map_into(LOs(exposed_sides2side.size(), 0), exposed_sides2side,
      side_class_ids_w, 1);
  if (classify_with & NODE_SETS) {
    int max_side_set_id = 0;
    if ((classify_with & SIDE_SETS) && side_set_ids.size()) {
      max_side_set_id =
          *std::max_element(side_set_ids.begin(), side_set_ids.end());
    }
    std::vector<int> node_set_ids(init_params.num_node_sets);
    CALL(ex_get_ids(file, EX_NODE_SET, node_set_ids.data()));
    for (size_t i = 0; i < node_set_ids.size(); ++i) {
      int nentries, ndist_factors;
      CALL(ex_get_set_param(
          file, EX_NODE_SET, node_set_ids[i], &nentries, &ndist_factors));
      if (verbose) {
        std::cout << "node set " << node_set_ids[i] << " has " << nentries
                  << " nodes\n";
      }
      if (ndist_factors) {
        std::cout << "Omega_h doesn't support distribution factors\n";
      }
      HostWrite<LO> h_set_nodes2nodes(nentries);
      CALL(ex_get_set(file, EX_NODE_SET, node_set_ids[i],
          h_set_nodes2nodes.data(), nullptr));
      auto set_nodes2nodes =
          subtract_from_each(LOs(h_set_nodes2nodes.write()), 1);
      auto nodes_are_in_set = mark_image(set_nodes2nodes, mesh->nverts());
      auto sides_are_in_set =
          mark_up_all(mesh, VERT, dim - 1, nodes_are_in_set);
      auto set_sides2side = collect_marked(sides_are_in_set);
      auto surface_id = node_set_ids[i] + max_side_set_id;
      if (verbose) {
        std::cout << "node set " << node_set_ids[i] << " will be surface "
                  << surface_id << '\n';
      }
      map_into(LOs(set_sides2side.size(), surface_id), set_sides2side,
          side_class_ids_w, 1);
      map_into(Read<I8>(set_sides2side.size(), I8(dim - 1)), set_sides2side,
          side_class_dims_w, 1);
    }
  }
  if (classify_with & SIDE_SETS) {
    for (size_t i = 0; i < side_set_ids.size(); ++i) {
      int nentries, ndist_factors;
      CALL(ex_get_set_param(
          file, EX_SIDE_SET, side_set_ids[i], &nentries, &ndist_factors));
      if (verbose) {
        std::cout << "side set " << side_set_ids[i] << " has " << nentries
                  << " sides\n";
      }
      if (ndist_factors && verbose) {
        std::cout << "Omega_h doesn't support distribution factors\n";
      }
      HostWrite<LO> h_set_sides2elem(nentries);
      HostWrite<LO> h_set_sides2local(nentries);
      CALL(ex_get_set(file, EX_SIDE_SET, side_set_ids[i],
          h_set_sides2elem.data(), h_set_sides2local.data()));
      auto set_sides2elem =
          subtract_from_each(LOs(h_set_sides2elem.write()), 1);
      auto set_sides2local = LOs(h_set_sides2local.write());
      auto elems2sides = mesh->ask_down(dim, dim - 1).ab2b;
      auto nsides_per_elem = simplex_degrees[dim][dim - 1];
      auto set_sides2side_w = Write<LO>(nentries);
      auto f2 = OMEGA_H_LAMBDA(LO set_side) {
        auto elem = set_sides2elem[set_side];
        auto local = side_exo2osh(dim, set_sides2local[set_side]);
        auto side = elems2sides[elem * nsides_per_elem + local];
        set_sides2side_w[set_side] = side;
      };
      parallel_for(nentries, f2);
      auto set_sides2side = LOs(set_sides2side_w);
      auto surface_id = side_set_ids[i];
      map_into(LOs(nentries, surface_id), set_sides2side, side_class_ids_w, 1);
      map_into(Read<I8>(nentries, I8(dim - 1)), set_sides2side,
          side_class_dims_w, 1);
    }
  }
  CALL(ex_close(file));
  auto elem_class_ids = LOs(elem_class_ids_w);
  auto side_class_ids = LOs(side_class_ids_w);
  auto side_class_dims = Read<I8>(side_class_dims_w);
  mesh->add_tag(dim, "class_id", 1, elem_class_ids);
  mesh->add_tag(dim - 1, "class_id", 1, side_class_ids);
  mesh->set_tag(dim - 1, "class_dim", side_class_dims);
  finalize_classification(mesh);
}

void write(
    std::string const& path, Mesh* mesh, bool verbose, int classify_with) {
  auto comp_ws = int(sizeof(Real));
  auto io_ws = comp_ws;
  auto mode = EX_CLOBBER | EX_MAPS_INT64_API;
  auto file = ex_create(path.c_str(), mode, &comp_ws, &io_ws);
  if (file < 0) Omega_h_fail("can't create Exodus file %s\n", path.c_str());
  auto title = "Omega_h " OMEGA_H_VERSION " Exodus Output";
  std::set<LO> region_set;
  auto dim = mesh->dim();
  auto elem_class_ids = mesh->get_array<LO>(dim, "class_id");
  auto h_elem_class_ids = HostRead<LO>(elem_class_ids);
  for (LO i = 0; i < h_elem_class_ids.size(); ++i) {
    region_set.insert(h_elem_class_ids[i]);
  }
  auto side_class_ids = mesh->get_array<LO>(dim - 1, "class_id");
  auto side_class_dims = mesh->get_array<I8>(dim - 1, "class_dim");
  auto h_side_class_ids = HostRead<LO>(side_class_ids);
  auto h_side_class_dims = HostRead<I8>(side_class_dims);
  std::set<LO> surface_set;
  for (LO i = 0; i < h_side_class_ids.size(); ++i) {
    if (h_side_class_dims[i] == I8(dim - 1)) {
      surface_set.insert(h_side_class_ids[i]);
    }
  }
  auto nelem_blocks = int(region_set.size());
  auto nside_sets =
      (classify_with | exodus::SIDE_SETS) ? int(surface_set.size()) : 0;
  auto nnode_sets =
      (classify_with | exodus::NODE_SETS) ? int(surface_set.size()) : 0;
  if (verbose) {
    std::cout << "init params for " << path << ":\n";
    std::cout << " Exodus ID " << file << '\n';
    std::cout << " comp_ws " << comp_ws << '\n';
    std::cout << " io_ws " << io_ws << '\n';
    std::cout << " Title " << title << '\n';
    std::cout << " num_dim " << dim << '\n';
    std::cout << " num_nodes " << mesh->nverts() << '\n';
    std::cout << " num_elem " << mesh->nelems() << '\n';
    std::cout << " num_elem_blk " << nelem_blocks << '\n';
    std::cout << " num_node_sets " << nnode_sets << '\n';
    std::cout << " num_side_sets " << nside_sets << '\n';
  }
  CALL(ex_put_init(file, title, dim, mesh->nverts(), mesh->nelems(),
      nelem_blocks, nnode_sets, nside_sets));
  Write<Real> coord_blk[3];
  for (Int i = 0; i < dim; ++i) coord_blk[i] = Write<Real>(mesh->nverts());
  auto coords = mesh->coords();
  auto f0 = OMEGA_H_LAMBDA(LO i) {
    for (Int j = 0; j < dim; ++j) coord_blk[j][i] = coords[i * dim + j];
  };
  parallel_for(mesh->nverts(), f0);
  HostRead<Real> h_coord_blk[3];
  for (Int i = 0; i < dim; ++i) h_coord_blk[i] = HostRead<Real>(coord_blk[i]);
  CALL(ex_put_coord(file, h_coord_blk[0].data(), h_coord_blk[1].data(),
      h_coord_blk[2].data()));
  auto all_conn = mesh->ask_elem_verts();
  for (auto block_id : region_set) {
    auto type_name = (dim == 3) ? "tetra4" : "tri3";
    auto elems_in_block = each_eq_to(elem_class_ids, block_id);
    auto block_elems2elem = collect_marked(elems_in_block);
    auto nblock_elems = block_elems2elem.size();
    if (verbose) {
      std::cout << "element block " << block_id << " has " << nblock_elems
                << " of type " << type_name << '\n';
    }
    CALL(ex_put_block(file, EX_ELEM_BLOCK, block_id, type_name, nblock_elems,
        dim + 1, 0, 0, 0));
    auto block_conn = unmap(block_elems2elem, all_conn, dim + 1);
    auto block_conn_ex = add_to_each(block_conn, 1);
    auto h_block_conn = HostRead<LO>(block_conn_ex);
    CALL(ex_put_conn(
        file, EX_ELEM_BLOCK, block_id, h_block_conn.data(), nullptr, nullptr));
  }
  if (classify_with) {
    for (auto set_id : surface_set) {
      auto sides_in_set = land_each(each_eq_to(side_class_ids, set_id),
          each_eq_to(side_class_dims, I8(dim - 1)));
      if (classify_with | exodus::SIDE_SETS) {
        auto set_sides2side = collect_marked(sides_in_set);
        auto nset_sides = set_sides2side.size();
        std::cout << "side set " << set_id << " has " << nset_sides
                  << " sides\n";
        auto sides2elems = mesh->ask_up(dim - 1, dim);
        Write<int> set_sides2elem(nset_sides);
        Write<int> set_sides2local(nset_sides);
        auto f1 = OMEGA_H_LAMBDA(LO set_side) {
          auto side = set_sides2side[set_side];
          auto se = sides2elems.a2ab[side];
          auto elem = sides2elems.ab2b[se];
          auto code = sides2elems.codes[se];
          auto which_down = code_which_down(code);
          set_sides2elem[set_side] = elem + 1;
          set_sides2local[set_side] = side_osh2exo(dim, which_down);
        };
        parallel_for(nset_sides, f1);
        auto h_set_sides2elem = HostRead<int>(set_sides2elem);
        auto h_set_sides2local = HostRead<int>(set_sides2local);
        CALL(ex_put_set_param(file, EX_SIDE_SET, set_id, nset_sides, 0));
        CALL(ex_put_set(file, EX_SIDE_SET, set_id, h_set_sides2elem.data(),
            h_set_sides2local.data()));
      }
      if (classify_with | exodus::NODE_SETS) {
        auto nodes_in_set = mark_down(mesh, dim - 1, VERT, sides_in_set);
        auto set_nodes2node = collect_marked(nodes_in_set);
        auto set_nodes2node_ex = add_to_each(set_nodes2node, 1);
        auto nset_nodes = set_nodes2node.size();
        std::cout << "node set " << set_id << " has " << nset_nodes
                  << " nodes\n";
        auto h_set_nodes2node = HostRead<LO>(set_nodes2node_ex);
        CALL(ex_put_set_param(file, EX_NODE_SET, set_id, nset_nodes, 0));
        CALL(ex_put_set(
            file, EX_NODE_SET, set_id, h_set_nodes2node.data(), nullptr));
      }
    }
  }
  CALL(ex_close(file));
}

#undef CALL

}  // end namespace exodus

}  // end namespace Omega_h
