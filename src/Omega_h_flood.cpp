#include "Omega_h_confined.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_loop.hpp"

namespace Omega_h {

Bytes mark_floodable_elements(Mesh* mesh) {
  auto edges_are_bridge = find_bridge_edges(mesh);
  auto vert_are_bridge = mark_down(mesh, EDGE, VERT, edges_are_bridge);
  auto elems_are_angle = find_angle_elems(mesh);
  auto verts_are_angle = mark_down(mesh, mesh->dim(), VERT, elems_are_angle);
  auto verts_can_flood = land_each(vert_are_bridge, verts_are_angle);
  auto elems_can_flood = mark_up(mesh, VERT, mesh->dim(), verts_can_flood);
  return elems_can_flood;
}

Bytes mark_flood_seeds(Mesh* mesh, AdaptOpts const& opts,
    Bytes elems_can_flood) {
  auto elem_quals = mesh->ask_qualities();
  /* I had thoughts of using much more complex criteria such as for pads
     comparing their distance to the minimum allowed edge length and for
     angles comparing the minimum possible quality given the forced dihedral
     angle to the minimum desired element quality.
     For now I'll just fall back on the lazy catch-all of low current quality.
     If misidentification of flood seeds becomes an issue, revisit the above
     ideas. */
  auto elems_are_lowqual = each_lt(elem_quals, opts.min_quality_desired);
  return land_each(elems_can_flood, elems_are_lowqual);
}

Bytes mark_seeded_flood_zones(Mesh* mesh, Bytes elems_can_flood, Bytes elems_are_seeded) {
  OMEGA_H_CHECK(elems_can_flood.size() == mesh->nelems());
  OMEGA_H_CHECK(marks.size() == mesh->nelems());
  auto dim = mesh->dim();
  auto e2s = mesh->ask_down(dim, dim - 1);
  auto s2e = mesh->ask_up(dim - 1, dim);
  auto side_class_dims = mesh->get_array<I8>(dim - 1, "class_dim");
  while (true) {
    auto elems_are_seeded_w = deep_copy(elems_are_seeded);
    auto f = OMEGA_H_LAMBDA(LO e) {
      if (!elems_can_flood[e]) return;
      for (Int ees = 0; ees < (dim + 1); ++ees) {
        auto s = e2s.ab2b[e * (dim + 1) + ees];
        if (side_class_dims[s] < dim) continue;
        auto se_begin = s2e.a2ab[s];
        auto se_end = s2e.a2ab[s + 1];
        for (LO se = se_begin; se < se_end; ++se) {
          auto oe = s2e.ab2b[se];
          if (oe == e) continue;
          if (!elems_can_flood[oe]) continue;
          if (marks[oe]) {
            elems_are_seeded_w[e] = Byte(1);
          }
        }
      }
    };
    parallel_for(mesh->nelems(), f);
    auto new_elems_are_seeded = Bytes(elems_are_seeded_w);
    new_elems_are_seeded = mesh->sync_array(dim, new_elems_are_seeded, 1);
    if (new_elems_are_seeded == elems_are_seeded) break;
    elems_are_seeded = new_elems_are_seeded;
  }
  return elems_are_seeded;
}

static
OMEGA_H_INLINE
void flood_update(I32& my_class_id, Real& my_density,
    I32 other_class_id, Real other_density) {
  if (other_class_id >= 0 &&
      (my_class_id == -1 ||
       other_density < my_density ||
       (other_density == my_density &&
        other_class_id < my_class_id))) {
    my_class_id = other_class_id;
    my_density = other_density;
  }
}

void flood_element_variables(Mesh* mesh, 
    Bytes elems_should_flood,
    Reals elem_densities,
    Read<I32>* p_elem_flood_class_ids,
    Reals* p_elem_flood_densities) {
  auto dim = mesh->dim();
  auto elem_class_ids = mesh->get_array<I32>(dim, "class_id");
  auto e2s = mesh->ask_down(dim, dim - 1);
  auto s2e = mesh->ask_up(dim - 1, dim);
  auto side_class_dims = mesh->get_array<I8>(dim - 1, "class_dim");
  // initially, no material is able to flood
  auto elem_flood_class_ids = Read<I32>(mesh->nelems(), I32(-1));
  auto elem_flood_densities = Reals(mesh->nelems(), ArithTraits<Reala>::max());
  while (true) {
    auto elem_flood_class_ids_w = deep_copy(elem_flood_class_ids);
    auto elem_flood_densities_w = deep_copy(elem_flood_densities);
    auto f = OMEGA_H_LAMBDA(LO e) {
      if (!elems_should_flood[e]) return;
      auto flood_class_id = elem_flood_class_ids[e];
      // once a flood material has been selected for an element,
      // it cannot be changed.
      // this will allow multiple flood fronts to advance in from
      // different boundaries and stop when they meet in the middle,
      // rather than forcing one material to fill the entire region.
      // certain whiteboard musings suggest that this will prevent
      // flip/flop livelock when two thin regions are trapped between objects.
      if (flood_class_id != -1) return;
      // however, we will accumulate all contributions across sides in
      // an associative way, to ensure the algorithm produces the
      // same answer regardless of element-local side ordering.
      auto flood_density = elem_flood_densities[e];
      for (Int ees = 0; ees < (dim + 1); ++ees) {
        auto s = e2s.ab2b[e * (dim + 1) + ees];
        // detect inter-material boundaries in the original model,
        // prior to flooding
        auto is_different_region = (side_class_dims[s] < dim);
        auto se_begin = s2e.a2ab[s];
        auto se_end = s2e.a2ab[s + 1];
        for (LO se = se_begin; se < se_end; ++se) {
          auto oe = s2e.ab2b[se];
          if (oe == e) continue;
          auto oe_should_flood = elems_should_flood[oe];
          // if the other element is across an old material boundary,
          // and the element on the other side is NOT part of another
          // flood region
          if (is_different_region && (!oe_should_flood)) {
            auto oe_class_id = elem_class_ids[oe];
            auto oe_density = elem_densities[oe];
            // we are allowed to start new flooding into the current flood
            // region from this other (original) material region
            flood_update(flood_class_id, flood_density,
                oe_class_id, oe_density);
          } else {
            auto oe_flood_class_id = elem_flood_class_ids[oe];
            auto oe_flood_density = elem_flood_densities[oe];
            // we are allowed to continue active flooding into the current
            // element
            flood_update(flood_class_id, flood_density,
                oe_flood_class_id, oe_flood_density);
          }
        }
      }
    };
    parallel_for(mesh->nelems(), f);
    auto new_elem_flood_class_ids = Read<I32>(elem_flood_class_ids_w);
    auto new_elem_flood_densities = Read<I32>(elem_flood_densities_w);
    new_elem_flood_class_ids = mesh->sync_array(dim, new_elem_flood_class_ids, 1);
    new_elem_flood_densities = mesh->sync_array(dim, new_elem_flood_densities, 1);
    if (new_elem_flood_class_ids == elem_flood_class_ids &&
        new_elem_flood_densities == elem_flood_densities) {
      break;
    }
    elem_flood_class_ids = new_elem_flood_class_ids;
    elem_flood_densities = new_elem_flood_densities;
  }
  *p_elem_flood_class_ids = elem_flood_class_ids;
  *p_elem_flood_densities = elem_flood_densities;
}
