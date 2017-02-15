#include <exodusII.h>
#include <algorithm>

#include "Omega_h.hpp"
#include "array.hpp"
#include "classify.hpp"
#include "internal.hpp"
#include "loop.hpp"
#include "map.hpp"
#include "mark.hpp"
#include "simplices.hpp"

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
        case 1: return 0;
        case 2: return 1;
        case 3: return 2;
      }
    case 3:
      switch (side) {
        case 1: return 1;
        case 2: return 2;
        case 3: return 3;
        case 4: return 0;
      }
  }
  return -1;
}

void read(std::string const& path, Mesh* mesh, bool verbose, int classify_with) {
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
  CALL(ex_get_coord(file,
        h_coord_blk[0].data(),
        h_coord_blk[1].data(),
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
    char elem_type[MAX_STR_LENGTH+1];
    int nentries;
    int nnodes_per_entry;
    int nedges_per_entry;
    int nfaces_per_entry;
    int nattr_per_entry;
    CALL(ex_get_block(file, EX_ELEM_BLOCK, block_ids[i], elem_type, &nentries,
          &nnodes_per_entry, &nedges_per_entry, &nfaces_per_entry,
          &nattr_per_entry));
    if (verbose) {
      std::cout << "block " << block_ids[i] << " has "
        << nentries << " elements of type " << elem_type << '\n';
    }
    if (!is_type_supported(dim, elem_type)) {
      Omega_h_fail("type %s is not supported for %dD !\n",
          elem_type, dim);
    }
    CHECK(nnodes_per_entry == dim + 1);
    std::vector<int> edge_conn(nentries * nedges_per_entry);
    std::vector<int> face_conn(nentries * nfaces_per_entry);
    CALL(ex_get_conn(file, EX_ELEM_BLOCK, block_ids[i],
          h_conn.data() + start, edge_conn.data(), face_conn.data()));
    auto region_id = block_ids[i];
    auto f0 = LAMBDA(LO entry) { elem_class_ids_w[start + entry] = region_id; };
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
  Write<I8> side_class_dims_w = deep_copy(mesh->get_array<I8>(dim - 1, "class_dim"));
  auto exposed_sides2side = collect_marked(sides_are_exposed);
  map_into(LOs(exposed_sides2side.size(), 0), exposed_sides2side,
      side_class_ids_w, 1);
  if (classify_with & NODE_SETS) {
    int max_side_set_id = 0;
    if ((classify_with & SIDE_SETS) && side_set_ids.size()) {
      max_side_set_id = *std::max_element(
          side_set_ids.begin(), side_set_ids.end());
    }
    std::vector<int> node_set_ids(init_params.num_node_sets);
    CALL(ex_get_ids(file, EX_NODE_SET, node_set_ids.data()));
    for (size_t i = 0; i < node_set_ids.size(); ++i) {
      int nentries, ndist_factors;
      CALL(ex_get_set_param(file, EX_NODE_SET, node_set_ids[i],
            &nentries, &ndist_factors));
      if (verbose) {
        std::cout << "node set " << node_set_ids[i] << " has "
          << nentries << " nodes\n";
      }
      if (ndist_factors) {
        std::cout << "Omega_h doesn't support distribution factors\n";
      }
      HostWrite<LO> h_set_nodes2nodes(nentries);
      CALL(ex_get_set(file, EX_NODE_SET, node_set_ids[i],
            h_set_nodes2nodes.data(), nullptr));
      auto set_nodes2nodes = subtract_from_each(LOs(h_set_nodes2nodes.write()), 1);
      auto nodes_are_in_set = mark_image(set_nodes2nodes, mesh->nverts());
      auto sides_are_in_set = mark_up_all(mesh, VERT, dim - 1, nodes_are_in_set);
      auto set_sides2side = collect_marked(sides_are_in_set);
      auto surface_id = node_set_ids[i] + max_side_set_id;
      if (verbose) {
        std::cout << "node set " << node_set_ids[i]
          << " will be surface " << surface_id << '\n';
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
      CALL(ex_get_set_param(file, EX_SIDE_SET, side_set_ids[i],
            &nentries, &ndist_factors));
      if (verbose) {
        std::cout << "side set " << side_set_ids[i] << " has "
          << nentries << " sides\n";
      }
      if (ndist_factors && verbose) {
        std::cout << "Omega_h doesn't support distribution factors\n";
      }
      HostWrite<LO> h_set_sides2elem(nentries);
      HostWrite<LO> h_set_sides2local(nentries);
      CALL(ex_get_set(file, EX_SIDE_SET, side_set_ids[i],
            h_set_sides2elem.data(), h_set_sides2local.data()));
      auto set_sides2elem = subtract_from_each(LOs(h_set_sides2elem.write()), 1);
      auto set_sides2local = LOs(h_set_sides2local.write());
      auto elems2sides = mesh->ask_down(dim, dim - 1).ab2b;
      auto nsides_per_elem = simplex_degrees[dim][dim - 1];
      auto set_sides2side_w = Write<LO>(nentries);
      auto f2 = LAMBDA(LO set_side) {
        auto elem = set_sides2elem[set_side];
        auto local = side_exo2osh(dim, set_sides2local[set_side]);
        auto side = elems2sides[elem * nsides_per_elem + local];
        set_sides2side_w[set_side] = side;
      };
      parallel_for(nentries, f2);
      auto set_sides2side = LOs(set_sides2side_w);
      auto surface_id = side_set_ids[i];
      map_into(LOs(nentries, surface_id), set_sides2side,
          side_class_ids_w, 1);
      map_into(Read<I8>(nentries, I8(dim - 1)), set_sides2side,
          side_class_dims_w, 1);
    }
  }
  CALL(ex_close(file));
  auto elem_class_ids = LOs(elem_class_ids_w);
  auto side_class_ids = LOs(side_class_ids_w);
  auto side_class_dims = Read<I8>(side_class_dims_w);
  mesh->add_tag(dim, "class_id", 1, OMEGA_H_INHERIT, OMEGA_H_DO_OUTPUT,
      elem_class_ids);
  mesh->add_tag(dim - 1, "class_id", 1, OMEGA_H_INHERIT, OMEGA_H_DO_OUTPUT,
      side_class_ids);
  mesh->set_tag(dim - 1, "class_dim", side_class_dims);
  finalize_classification(mesh);
}

#undef CALL

} // end namespace exodus

} // end namespace Omega_h
