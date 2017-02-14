#include <exodusII.h>

#include "Omega_h.hpp"
#include "internal.hpp"

namespace Omega_h {

#define CALL(f) CHECK((f) >= 0)

namespace exodus {

static char const* const elem_types[4] = {
  "",
  "",
  "tri3",
  "tetra4"
};

void read(std::string const& path, Mesh* mesh, bool verbose) {
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
  auto dim = init_params.num_dim;
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
    CHECK(std::string(elem_type) == elem_types[dim]);
    CHECK(nnodes_per_entry == dim + 1);
    std::vector<int> edge_conn(nentries * nedges_per_entry);
    std::vector<int> face_conn(nentries * nfaces_per_entry);
    CALL(ex_get_conn(file, EX_ELEM_BLOCK, block_ids[i],
          h_conn.data() + start, edge_conn.data(), face_conn.data()));
    start += nentries * nnodes_per_entry;
  }
  CHECK(start == init_params.num_elem * (dim + 1));
  auto conn = LOs(h_conn.write());
  build_from_elems_and_coords(mesh, dim, conn, coords);
//if (init_params.num_node_maps) {
//  std::vector<int> node_map_ids(init_params.num_node_maps);
//  CALL(ex_get_ids(file, EX_NODE_MAP, node_map_ids.data()));
//}
}

#undef CALL

} // end namespace exodus

} // end namespace Omega_h
