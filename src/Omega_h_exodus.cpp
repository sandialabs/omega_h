#include "Omega_h.hpp"

#include <exodusII.h>

namespace Omega_h {

#define CALL(f) CHECK((f) >= 0)

namespace exodus {

static char const* const elem_types[4] = {
  "",
  "",
  "TRI3",
  "TETRA4"
};

void read(std::string const& path, Mesh* mesh, bool verbose) {
  auto comp_ws = int(sizeof(Real));
  int io_ws;
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
  HostRead<Real> h_coord_blk[3];
  for (Int i = 0; i < init_params.num_dim; ++i) {
    h_coord_blk[i] = HostRead<Real>(init_params.num_nodes);
  }
  CALL(ex_get_coord(file,
        h_coord_blk[0].data(),
        h_coord_blk[1].data(),
        h_coord_blk[2].data()));
  std::vector<int> block_ids(init_params.num_elem_blk);
  std::vector<int> block_conns(init_params.num_elem_blk);
  CALL(ex_get_ids(file, EX_ELEM_BLOCK, block_ids.data()));
  for (Int i = 0; i < block_ids.size(); ++i) {
    char elem_type[MAX_STR_LENGTH+1];
    int nentries;
    int nnodes_per_entry;
    int nedges_per_entry;
    int nfaces_per_entry;
    int nattr_per_entry;
    CALL(ex_get_block(file, EX_ELEM_BLOCK, block_ids[i],
          elem_type, &nnodes_per_entry, &nedges_per_entry,
          &nfaces_per_entry, &nattr_per_entry));
    CHECK(std::string(elem_type) == elem_types[init_params.num_dim]);
    CHECK(nnodes_per_entry == init_params.num_dim + 1);
    block_conns[i].resize(nentries * nnodes_per_entry);
    std::vector<int> edge_conn(nentries * nedges_per_entry);
    std::vector<int> face_conn(nentries * nfaces_per_entry);
    CALL(ex_get_conn(file, EX_ELEM_BLOCK, block_ids[i],
          block_conns[i].data(), edge_conn.data(), face_conn.data()));
  }
}

void write(std::string const& path, Mesh* mesh) {
}

#undef CALL

} // end namespace Omega_h
