#ifndef OMEGA_H_STK_HPP
#define OMEGA_H_STK_HPP

namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class BulkData; } }

namespace Omega_h {

  namespace stk_mesh {
    void read_from_stk(Mesh* mesh,
                       const stk::mesh::BulkData & stk_mesh,
                       const std::vector<std::pair<std::string,std::set<std::string>>> & merging_block_names);
    void read_from_stk(Mesh* mesh,
                       const stk::mesh::BulkData & stk_mesh,
                       const std::vector<stk::mesh::Entity> & stk_elems,
                       const std::vector<std::pair<std::string,std::set<std::string>>> & merging_block_names);
  }  // namespace stk_mesh

}  // namespace Omega_h

#endif
