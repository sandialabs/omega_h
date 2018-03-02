#include <Omega_h_array_ops.hpp>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_mesh.hpp>

namespace Omega_h {

LOs count_amr_new_ents(Mesh* mesh, Read<Byte> elem_mark) {
  Write<LO> num_new_ents(4, 0);
  auto dim = mesh->dim();
  for (Int i = 1; i <= dim; ++i) {
    auto dim_mark = mark_down(mesh, dim, i, elem_mark);
    for (Int j = 0; j <= i; ++j) {
      auto deg = hypercube_split_degree(i, j);
      auto nsplit = get_sum(dim_mark);
      num_new_ents[j] += deg * nsplit;
    }
  }
  return num_new_ents;
}

}
