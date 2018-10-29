#include <Omega_h_amr_topology.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_mesh.hpp>

namespace Omega_h {

namespace amr {

static Bytes mark_leaf_down(Mesh* mesh, Int mod_dim, Bytes elems_are_marked) {
  auto elem_dim = mesh->dim();
  auto is_mod_dim_leaf = mesh->ask_leaves(mod_dim);
  auto dim_mark = mark_down(mesh, elem_dim, mod_dim, elems_are_marked);
  return land_each(is_mod_dim_leaf, dim_mark);
}

void mark_refined(Mesh* mesh, Bytes elems_are_marked) {
  auto elem_dim = mesh->dim();
  for (Int mod_dim = 0; mod_dim <= elem_dim; ++mod_dim) {
    auto mark = mark_leaf_down(mesh, mod_dim, elems_are_marked);
    mesh->add_tag<Omega_h::Byte>(mod_dim, "refine", 1, mark);
  }
}

Few<LO, 4> count_refined(Mesh* mesh) {
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

LOs get_refined_topology(Mesh* mesh, Int child_dim, LO num_children,
    Few<LOs, 4> mods2mds, Few<LOs, 4> mds2mods, Few<LOs, 4> mods2midverts,
    LOs old_verts2new_verts) {
  Int spatial_dim = mesh->dim();
  Int num_verts_per_child = hypercube_degree(child_dim, 0);
  Write<LO> child_verts(num_children * num_verts_per_child);
  LO offset = 0;
  Few<Children, 4> old_parents2child_verts;
  for (Int lowd = 1; lowd <= spatial_dim; ++lowd) {
    old_parents2child_verts[lowd] = mesh->ask_children(lowd, 0);
  }
  for (Int d = child_dim; d <= spatial_dim; ++d) {
    Few<LOs, 4> mds2lows;
    for (Int lowd = 0; lowd <= d; ++lowd) {
      mds2lows[lowd] = mesh->ask_graph(d, lowd).ab2b;
    }
    LO num_mods = mods2midverts[d].size();
    Int num_child_per_mod = hypercube_split_degree(d, child_dim);
    auto mod_loop = OMEGA_H_LAMBDA(LO mod) {
      LO md = mods2mds[d][mod];
      for (Int child = 0; child < num_child_per_mod; ++child) {
        for (Int vert = 0; vert < num_verts_per_child; ++vert) {
          auto low = hypercube_split_template(d, child_dim, child, vert);
          Int num_lows_per_mod = hypercube_degree(d, low.dim);
          LO low_id = mds2lows[low.dim][md * num_lows_per_mod + low.which_down];
          LO midvert;
          LO low_adj_mod = mds2mods[low.dim][low_id];
          if (low_adj_mod == -1) {
            auto begin = old_parents2child_verts[low.dim].a2ab[low_id];
            auto end = old_parents2child_verts[low.dim].a2ab[low_id + 1];
            OMEGA_H_CHECK((end - begin) == 1);
            auto old_midvert = old_parents2child_verts[low.dim].ab2b[begin];
            midvert = old_verts2new_verts[old_midvert];
          } else {
            midvert = mods2midverts[low.dim][low_adj_mod];
          }
          LO idx = offset +
                   (mod * num_child_per_mod + child) * num_verts_per_child +
                   vert;
          child_verts[idx] = midvert;
        }
      }
    };
    parallel_for(num_mods, mod_loop, "get_amr_topology");
    offset += num_mods * num_child_per_mod * num_verts_per_child;
  }
  return child_verts;
}

}  // namespace amr

}  // namespace Omega_h
