#include <Omega_h_array_ops.hpp>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_mesh.hpp>

namespace Omega_h {

void mark_amr(Mesh* mesh, Read<Byte> elem_mark) {
  auto dim = mesh->dim();
  for (Int i = 1; i <= dim; ++i) {
    auto dim_mark = mark_down(mesh, dim, i, elem_mark);
    mesh->add_tag<Omega_h::Byte>(i, "refine", 1, dim_mark);
  }
}

Few<Real, 4> count_amr(Mesh* mesh) {
  auto dim = mesh->dim();
  Few<Real, 4> num_ents({0,0,0,0});
  for (Int i = 1; i <=dim; ++i) {
    auto dim_tag = mesh->get_tag<Byte>(i, "refine");
    for (Int j = 0; j <= i; ++j) {
      auto deg = hypercube_split_degree(i, j);
      auto nsplit = get_sum(dim_tag->array());
      num_ents[j] += deg * nsplit;
    }
  }
  return num_ents;
}

}
