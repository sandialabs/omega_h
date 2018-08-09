#include "Omega_h_flood.hpp"

#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_confined.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"

#include <iostream>
#include "Omega_h_file.hpp"

namespace Omega_h {

Bytes mark_flood_zones(Mesh* mesh, FloodOpts const&) {
  auto dim = mesh->dim();
  auto edges_are_bridge = find_bridge_edges(mesh);
  auto elems_adj_bridge = mark_adj(mesh, EDGE, mesh->dim(), edges_are_bridge);
  auto elems_can_seed = elems_adj_bridge;
  auto elem_qualities = mesh->ask_qualities();
  /* this is a totally useless standin to just compile while we work on the new
   * stuff */
  auto elems_low_qual = each_lt(elem_qualities, 0.5);
  auto elems_are_seeded = land_each(elems_low_qual, elems_can_seed);
  auto verts_are_seeds = mark_adj(mesh, dim, VERT, elems_are_seeded);
  auto elems_can_flood = mark_adj(mesh, VERT, dim, verts_are_seeds);
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
          if (elems_are_seeded[oe]) {
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

static OMEGA_H_INLINE void flood_update(I32& my_class_id, Real& my_density,
    I32 other_class_id, Real other_density) {
  if (other_class_id >= 0 &&
      (my_class_id == -1 || other_density < my_density ||
          (other_density == my_density && other_class_id < my_class_id))) {
    my_class_id = other_class_id;
    my_density = other_density;
  }
}

void flood_element_variables(Mesh* mesh, Bytes elems_should_flood,
    Reals elem_densities, Read<I32>* p_elem_flood_class_ids,
    Reals* p_elem_flood_densities) {
  auto dim = mesh->dim();
  OMEGA_H_CHECK(mesh->owners_have_all_upward(dim - 1));
  auto elem_class_ids = mesh->get_array<ClassId>(dim, "class_id");
  auto e2s = mesh->ask_down(dim, dim - 1);
  auto s2e = mesh->ask_up(dim - 1, dim);
  auto side_class_dims = mesh->get_array<I8>(dim - 1, "class_dim");
  // initially, no material is able to flood
  auto elem_flood_class_ids = Read<I32>(mesh->nelems(), I32(-1));
  auto elem_flood_densities = Reals(mesh->nelems(), ArithTraits<Real>::max());
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
            flood_update(
                flood_class_id, flood_density, oe_class_id, oe_density);
          } else {
            auto oe_flood_class_id = elem_flood_class_ids[oe];
            auto oe_flood_density = elem_flood_densities[oe];
            // we are allowed to continue active flooding into the current
            // element
            flood_update(flood_class_id, flood_density, oe_flood_class_id,
                oe_flood_density);
          }
        }
      }
      elem_flood_class_ids_w[e] = flood_class_id;
      elem_flood_densities_w[e] = flood_density;
    };
    parallel_for(mesh->nelems(), f);
    auto new_elem_flood_class_ids = Read<I32>(elem_flood_class_ids_w);
    auto new_elem_flood_densities = Reals(elem_flood_densities_w);
    new_elem_flood_class_ids =
        mesh->sync_array(dim, new_elem_flood_class_ids, 1);
    new_elem_flood_densities =
        mesh->sync_array(dim, new_elem_flood_densities, 1);
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

// tries to extend existing boundary tags onto
// the newly formed interfaces created by flooding
void flood_class_ids(Mesh* mesh, Int ent_dim) {
  auto class_ids = mesh->get_array<ClassId>(ent_dim, "class_id");
  auto class_dims = mesh->get_array<I8>(ent_dim, "class_dim");
  auto low_class_dims = mesh->get_array<I8>(ent_dim - 1, "class_dim");
  while (true) {
    auto class_ids_w = deep_copy(class_ids);
    auto h2l = mesh->ask_down(ent_dim, ent_dim - 1);
    auto l2h = mesh->ask_up(ent_dim - 1, ent_dim);
    auto f = OMEGA_H_LAMBDA(LO h) {
      auto class_dim = class_dims[h];
      if (class_dim != ent_dim) return;
      auto class_id = class_ids[h];
      if (class_id != -1) return;
      for (LO hl = h * (ent_dim + 1); hl < (h + 1) * (ent_dim + 1); ++hl) {
        auto l = h2l.ab2b[hl];
        auto l_class_dim = low_class_dims[l];
        if (l_class_dim < ent_dim) continue;
        for (LO lh = l2h.a2ab[l]; lh < l2h.a2ab[l + 1]; ++lh) {
          auto oh = l2h.ab2b[lh];
          auto oh_class_dim = class_dims[oh];
          if (oh_class_dim != class_dim) continue;
          auto oh_class_id = class_ids[oh];
          if (oh_class_id == -1) continue;
          if (class_id == -1 || oh_class_id < class_id) {
            class_id = oh_class_id;
          }
        }
      }
      class_ids_w[h] = class_id;
    };
    parallel_for(mesh->nents(ent_dim), f);
    auto new_class_ids = Read<ClassId>(class_ids_w);
    new_class_ids = mesh->sync_array(ent_dim, new_class_ids, 1);
    if (new_class_ids == class_ids) {
      break;
    }
    class_ids = new_class_ids;
  }
  mesh->set_tag(ent_dim, "class_id", class_ids);
}

void flood_classification(Mesh* mesh, Bytes elems_did_flood) {
  auto dim = mesh->dim();
  for (Int ent_dim = VERT; ent_dim < dim; ++ent_dim) {
    auto ents_in_flood_closure = mark_down(mesh, dim, ent_dim, elems_did_flood);
    std::cout << "clearing classification for "
              << get_sum(ents_in_flood_closure) << " entities of dimension "
              << ent_dim << '\n';
    auto class_ids_w = deep_copy(mesh->get_array<ClassId>(ent_dim, "class_id"));
    auto class_dims_w = deep_copy(mesh->get_array<I8>(ent_dim, "class_dim"));
    auto f = OMEGA_H_LAMBDA(LO i) {
      if (ents_in_flood_closure[i]) {
        class_ids_w[i] = -1;
        class_dims_w[i] = I8(dim);
      }
    };
    parallel_for(mesh->nents(ent_dim), f);
    mesh->set_tag(ent_dim, "class_id", Read<ClassId>(class_ids_w));
    mesh->set_tag(ent_dim, "class_dim", Read<I8>(class_dims_w));
  }
  for (Int ent_dim = dim - 1; ent_dim >= VERT; --ent_dim) {
    bool relaxed = true;
    project_classification(mesh, ent_dim, relaxed);
    if (ent_dim > VERT) {
      // we have to project one further down,
      // because flood_class_ids() looks at those
      // entities to decide how far IDs can flood
      project_classification(mesh, ent_dim - 1, relaxed);
      flood_class_ids(mesh, ent_dim);
    }
  }
}

void flood(Mesh* mesh, FloodOpts const& opts) {
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto dim = mesh->dim();
  auto elems_should_flood = mark_flood_zones(mesh, opts);
  auto elem_densities = mesh->get_array<Real>(dim, opts.density_name);
  Read<I32> elem_flood_class_ids;
  Reals elem_flood_densities;
  flood_element_variables(mesh, elems_should_flood, elem_densities,
      &elem_flood_class_ids, &elem_flood_densities);
  auto elem_class_ids = mesh->get_array<ClassId>(dim, "class_id");
  auto elem_class_ids_w = Write<ClassId>(mesh->nelems());
  auto elems_did_flood_w = Write<I8>(mesh->nelems());
  auto elem_densities_w = Write<Real>(mesh->nelems());
  auto f = OMEGA_H_LAMBDA(LO e) {
    if (elem_flood_class_ids[e] == -1) {
      elems_did_flood_w[e] = 0;
      elem_class_ids_w[e] = elem_class_ids[e];
      elem_densities_w[e] = elem_densities[e];
    } else {
      elems_did_flood_w[e] = 1;
      elem_class_ids_w[e] = elem_flood_class_ids[e];
      elem_densities_w[e] = elem_flood_densities[e];
    }
  };
  parallel_for(mesh->nelems(), f);
  auto elems_did_flood = Bytes(elems_did_flood_w);
  std::cout << get_sum(elems_did_flood) << " elements flooded\n";
  mesh->set_tag(dim, "class_id", Read<ClassId>(elem_class_ids_w));
  mesh->set_tag(dim, opts.density_name, Reals(elem_densities_w));
  flood_classification(mesh, elems_did_flood);
}

}  // end namespace Omega_h
