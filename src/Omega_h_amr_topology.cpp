#include <Omega_h_amr_topology.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_loop.hpp>
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

Few<LO, 4> count_amr(Mesh* mesh) {
  auto dim = mesh->dim();
  Few<LO, 4> num_ents({0, 0, 0, 0});
  for (Int i = 1; i <= dim; ++i) {
    auto dim_mark = mesh->get_array<Byte>(i, "refine");
    for (Int j = 0; j <= i; ++j) {
      auto deg = hypercube_split_degree(i, j);
      auto nsplit = get_sum(dim_mark);
      num_ents[j] += deg * nsplit;
    }
  }
  return num_ents;
}

LOs get_amr_topology(Mesh* mesh, Int child_dim, LO num_children,
    Few<LOs, 4> parents2mds, Few<LOs, 4> mds2parents,
    Few<LOs, 4> parents2midverts) {
  Int spatial_dim = mesh->dim();
  Int num_verts_per_child = hypercube_degree(child_dim, 0);
  Write<LO> child_verts(num_children * num_verts_per_child);
  for (Int dim = 0; dim <= spatial_dim; ++dim) {
    Int offset = 0;
    Few <LOs, 4> mds2lows;
    for (Int low_dim = 0; low_dim <= 3; ++low_dim)
      mds2lows[low_dim] = mesh->ask_down(dim, low_dim).ab2b;
    LO num_parents = parents2midverts[dim].size();
    Int num_child_per_parent = hypercube_split_degree(dim, child_dim);
    auto parent_loop = OMEGA_H_LAMBDA(LO parent) {
      LO md = parents2mds[dim][parent];
      for (Int child = 0; child < num_child_per_parent; ++child) {
        for (Int vert = 0; vert < num_verts_per_child; ++vert) {
          auto low = hypercube_split_template(dim, child_dim, child, vert);
          Int num_lows_per_parent = hypercube_degree(dim, low.dim);
          LO low_gid = mds2lows[low.dim][md * num_lows_per_parent + low.which_down];
          LO low_adj_parent = mds2parents[dim][low_gid];
          LO midvert = parents2midverts[dim][low_adj_parent];
          LO idx = offset + (parent * num_child_per_parent + child) * num_verts_per_child + vert;
          child_verts[idx] = midvert;
        }
      }
    };
    parallel_for(num_parents, parent_loop);
    offset += num_parents * num_child_per_parent * num_verts_per_child;
  }
  return child_verts;
}

}  // namespace Omega_h
