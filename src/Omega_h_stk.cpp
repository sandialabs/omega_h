#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <algorithm>
#include <set>
#include <map>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_stk.hpp"

namespace Omega_h {


namespace stk_mesh {

static int get_element_part_id(const stk::mesh::BulkData& mesh, stk::mesh::Entity elem)
{
  OMEGA_H_CHECK(mesh.entity_rank(elem) == stk::topology::ELEMENT_RANK);
  const stk::mesh::Part * elem_io_part = nullptr;

  for(auto && part : mesh.bucket(elem).supersets())
  {
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK && part->subsets().empty() && part->topology() != stk::topology::INVALID_TOPOLOGY)
    {
      // there should only be one element rank part without subsets with topology on the element
      OMEGA_H_CHECK(nullptr == elem_io_part);
      elem_io_part = part;
    }
  }
  OMEGA_H_CHECK(NULL != elem_io_part);
  return elem_io_part->id();
}

static std::vector<stk::mesh::Entity> get_sorted_elem_nodes(const stk::mesh::BulkData & stk_mesh, const std::vector<stk::mesh::Entity> & stk_elems)
{
  stk::mesh::EntityLess less(stk_mesh);
  auto dim = int(stk_mesh.mesh_meta_data().spatial_dimension());
  const int nodes_per_elem = (dim == 3) ? 4 : 3;

  std::set<stk::mesh::Entity, stk::mesh::EntityLess> nodes(less);
  for (auto && elem : stk_elems)
  {
    const stk::mesh::Entity * elem_nodes = stk_mesh.begin_nodes(elem);
    OMEGA_H_CHECK((int)stk_mesh.num_nodes(elem) == nodes_per_elem);
    for (int i=0; i<nodes_per_elem; ++i)
    {
      nodes.insert(elem_nodes[i]);
    }
  }

  std::vector<stk::mesh::Entity> elem_nodes(nodes.begin(), nodes.end());
  return elem_nodes;
}

static LO get_node_index(const stk::mesh::BulkData & stk_mesh, stk::mesh::Entity node, const std::vector<stk::mesh::Entity> & sorted_nodes)
{
  auto elem_node = std::lower_bound(sorted_nodes.begin(), sorted_nodes.end(), node, stk::mesh::EntityLess(stk_mesh));
  OMEGA_H_CHECK(elem_node != sorted_nodes.end() && *elem_node == node);
  return LO(std::distance(sorted_nodes.begin(), elem_node));
}

static LOs get_elem_node_indices(const stk::mesh::BulkData & stk_mesh, const std::vector<stk::mesh::Entity> & stk_elems, const std::vector<stk::mesh::Entity> & sorted_nodes)
{
  stk::mesh::EntityLess less(stk_mesh);
  auto dim = int(stk_mesh.mesh_meta_data().spatial_dimension());
  const int nodes_per_elem = (dim == 3) ? 4 : 3;

  HostWrite<LO> elem_node_indices(LO(stk_elems.size() * nodes_per_elem));

  LO index = 0;
  for (auto && elem : stk_elems)
  {
    const stk::mesh::Entity * elem_nodes = stk_mesh.begin_nodes(elem);
    OMEGA_H_CHECK((int)stk_mesh.num_nodes(elem) == nodes_per_elem);
    for (int i=0; i<nodes_per_elem; ++i)
    {
      elem_node_indices[index++] = get_node_index(stk_mesh, elem_nodes[i], sorted_nodes);
    }
  }

  return elem_node_indices.write();
}

static LOs get_elem_block_ids_compress(const stk::mesh::BulkData & stk_mesh,
                                       const std::vector<stk::mesh::Entity> & stk_elems,
                                       const std::vector<std::set<int>> & merging_blocks_original_IDs,
                                       std::map<LO,LO>& unique_ids_map)
{
  HostWrite<LO> elem_block_ids(LO(stk_elems.size()));
  LO unique_count = 0;

  // data for merging blocks
  std::vector<LO> new_merged_unique_ID(merging_blocks_original_IDs.size());
  std::vector<bool> block_set_processed(merging_blocks_original_IDs.size(), false);

  LO index = 0;
  for (auto && elem : stk_elems)
  {
    // check if this element's block needs to be added to unique_ids_map
    auto const original_id = get_element_part_id(stk_mesh, elem);
    bool add_unique = true;
    for (auto const& a_pair : unique_ids_map) {
      if (a_pair.first == original_id) {
        add_unique = false;
        break;
      }
    }

    if (add_unique) {

      // check if element is in a set of blocks to be merged
      size_t block_set_index = 0;
      bool block_is_merging = false;
      for (size_t b=0; b<merging_blocks_original_IDs.size(); ++b) {
        std::set<int>::const_iterator iter = merging_blocks_original_IDs[b].find(original_id);
        if (iter != merging_blocks_original_IDs[b].end()) {
          if (!block_set_processed[b]) {
            new_merged_unique_ID[b] = unique_count++;
            block_set_processed[b] = true;
          }
          block_set_index = b;
          block_is_merging = true;
          break;
        }
      }

      if (block_is_merging)
        unique_ids_map.insert({original_id,new_merged_unique_ID[block_set_index]});
      else
        unique_ids_map.insert({original_id,unique_count++});
    }
    elem_block_ids[index++] = unique_ids_map[original_id];
  }
  return elem_block_ids.write();
}

static LOs get_elem_block_ids(const stk::mesh::BulkData & stk_mesh, const std::vector<stk::mesh::Entity> & stk_elems)
{
  HostWrite<LO> elem_block_ids(LO(stk_elems.size()));

  LO index = 0;
  for (auto && elem : stk_elems)
  {
    elem_block_ids[index++] = get_element_part_id(stk_mesh, elem);
  }

  return elem_block_ids.write();
}

static Reals get_node_coords(const stk::mesh::BulkData & stk_mesh, const std::vector<stk::mesh::Entity> & nodes)
{
  auto dim = int(stk_mesh.mesh_meta_data().spatial_dimension());
  const stk::mesh::FieldBase * coords_field(stk_mesh.mesh_meta_data().coordinate_field());

  HostWrite<Real> node_coords(nodes.size()*dim);
  int index = 0;
  for (auto && node : nodes)
  {
    const double * coord_data = static_cast<double *>(stk::mesh::field_data(*coords_field, node));
    for (Int j = 0; j < dim; ++j) {
      node_coords[index++] = coord_data[j];
    }
  }

  return node_coords.write();
}

static int get_side_id(const stk::mesh::BulkData & stk_mesh, stk::mesh::Entity side)
{
  int best_id = -1;
  if (stk_mesh.is_valid(side))
  {
    for(auto && part : stk_mesh.bucket(side).supersets())
    {
      if (part->primary_entity_rank() == stk::topology::FACE_RANK && part->id() > 0)
        return int(part->id());
    }
  }
  return best_id;
}

static std::vector<stk::mesh::Entity> get_side_entities(Mesh* mesh, const stk::mesh::BulkData & stk_mesh, const std::vector<stk::mesh::Entity> & sorted_nodes)
{
  auto dim = int(stk_mesh.mesh_meta_data().spatial_dimension());
  auto side_verts = mesh->ask_verts_of(dim - 1);
  HostRead<LO> host_side_verts = HostRead<LO>(side_verts);
  auto num_side_verts = element_degree(mesh->family(), dim - 1, VERT);

  std::vector<stk::mesh::Entity> side_vec;
  std::vector<stk::mesh::Entity> side_nodes(num_side_verts);

  std::vector<stk::mesh::Entity> side_entities(mesh->nents(dim - 1));
  for (size_t i=0; i<side_entities.size(); ++i)
  {
    for (int v=0; v<num_side_verts; ++v)
    {
      side_nodes[v] = sorted_nodes[host_side_verts[i*num_side_verts+v]];
    }
    stk::mesh::get_entities_through_relations(stk_mesh, side_nodes, stk_mesh.mesh_meta_data().side_rank(), side_vec);
    OMEGA_H_CHECK(side_vec.size() <= 1);
    if (!side_vec.empty()) side_entities[i] = side_vec[0];
  }
  return side_entities;
}

static void classify_elements(Mesh * mesh,
                              const stk::mesh::BulkData & stk_mesh,
                              const std::vector<stk::mesh::Entity> & stk_elems,
                              const std::vector<std::pair<std::string,std::set<std::string>>> & merging_block_names)
{
  auto dim = int(stk_mesh.mesh_meta_data().spatial_dimension());
  classify_elements(mesh);

  // store ID's of blocks to be merged by block set
  std::vector<stk::mesh::Part*> const& myparts = stk_mesh.mesh_meta_data().get_parts();
  std::vector<std::set<int>> merging_blocks_original_IDs;
  std::map<std::string,std::string> old2new_names;

  if (merging_block_names.size() > 0) {
    std::set<int> merge_IDs;
    for(size_t b=0; b<merging_block_names.size(); ++b) {
      merge_IDs.clear();
      for (size_t k=0; k<myparts.size(); ++k) {
        if (myparts[k]->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
            myparts[k]->id() > 0) {

          std::set<std::string>::const_iterator iter = merging_block_names[b].second.find(myparts[k]->name());
          if (iter != merging_block_names[b].second.end()) {
            merge_IDs.insert(myparts[k]->id());
            old2new_names[myparts[k]->name()] = merging_block_names[b].first;
          }
        }
      }
      merging_blocks_original_IDs.push_back(merge_IDs);
    }
  }

  std::map<LO,LO> unique_ids_map;
  LOs elem_class_ids = get_elem_block_ids_compress(stk_mesh, stk_elems, merging_blocks_original_IDs, unique_ids_map);
  mesh->add_tag(dim, "class_id", 1, elem_class_ids);
  std::set<std::string> new_names_added;
  for (size_t k=0; k<myparts.size(); ++k) {

    bool id_used_by_elements = false;
    std::map<LO,LO>::const_iterator iter = unique_ids_map.find(myparts[k]->id());
    if (iter != unique_ids_map.end())
      id_used_by_elements = true;
    if (myparts[k]->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
        myparts[k]->id() > 0 &&
        id_used_by_elements) {

      // process blocks with changing names due to merging
      bool name_is_changing = false;
      if (old2new_names.size() > 0) {
        std::map<std::string,std::string>::const_iterator old2new_names_iter = old2new_names.find(myparts[k]->name());
        if (old2new_names_iter != old2new_names.end()) {
          name_is_changing = true;
          std::set<std::string>::iterator new_names_iter = new_names_added.find(old2new_names_iter->second);
          if (new_names_iter == new_names_added.end()) { // not found so add it
            int the_id = int(myparts[k]->id());
            int the_new_id = unique_ids_map[the_id];
            mesh->class_sets[old2new_names_iter->second].push_back({dim, the_new_id});
            new_names_added.insert(old2new_names_iter->second);
          }
        }
      }

      // process blocks whose names are not changing
      if (!name_is_changing) {
        int the_id = int(myparts[k]->id());
        int the_new_id = unique_ids_map[the_id];
        mesh->class_sets[myparts[k]->name()].push_back({dim, the_new_id});
      }
    }
  }
}

static void classify_sides(Mesh * mesh, const stk::mesh::BulkData & stk_mesh, const std::vector<stk::mesh::Entity> & sorted_nodes)
{
  auto dim = int(stk_mesh.mesh_meta_data().spatial_dimension());
  LO num_sides = mesh->nents(dim - 1);
  Write<LO> side_class_ids_w(num_sides, -1);
  auto sides_are_exposed = mark_exposed_sides(mesh);
  classify_sides_by_exposure(mesh, sides_are_exposed);
  Write<I8> side_class_dims_w =
      deep_copy(mesh->get_array<I8>(dim - 1, "class_dim"));
  auto exposed_sides2side = collect_marked(sides_are_exposed);
  map_into(LOs(exposed_sides2side.size(), 0), exposed_sides2side,
      side_class_ids_w, 1);
  auto host_side_class_dims_w = HostWrite<I8>(side_class_dims_w);
  auto host_side_class_ids_w = HostWrite<LO>(side_class_ids_w);

  const std::vector<stk::mesh::Entity> side_entities = get_side_entities(mesh, stk_mesh, sorted_nodes);

  for (LO i=0; i<num_sides; ++i)
  {
    int side_id = get_side_id(stk_mesh, side_entities[i]);
    if (side_id >= 0)
    {
      host_side_class_ids_w[i] = side_id;
      host_side_class_dims_w[i] = I8(dim - 1);
    }
  }

  auto side_class_ids = LOs(host_side_class_ids_w);
  auto side_class_dims = Read<I8>(host_side_class_dims_w);
  mesh->add_tag(dim - 1, "class_id", 1, side_class_ids);
  mesh->set_tag(dim - 1, "class_dim", side_class_dims);

  std::vector<stk::mesh::Part*> const& myparts = stk_mesh.mesh_meta_data().get_parts();
  for (size_t k=0; k<myparts.size(); ++k)
    if (myparts[k]->primary_entity_rank() == stk::topology::FACE_RANK &&
        myparts[k]->id() > 0)
      mesh->class_sets[myparts[k]->name()].push_back({dim-1, int(myparts[k]->id())});
}

void read_from_stk(Mesh* mesh,
                   const stk::mesh::BulkData & stk_mesh,
                   const std::vector<stk::mesh::Entity> & stk_elems,
                   const std::vector<std::pair<std::string,std::set<std::string>>> & merging_block_names)
{
  begin_code("stk::clone_from_stk");

  const std::vector<stk::mesh::Entity> sorted_nodes = get_sorted_elem_nodes(stk_mesh, stk_elems);

  LOs conn = get_elem_node_indices(stk_mesh, stk_elems, sorted_nodes);
  Reals coords = get_node_coords(stk_mesh, sorted_nodes);

  auto dim = int(stk_mesh.mesh_meta_data().spatial_dimension());
  build_from_elems_and_coords(mesh, OMEGA_H_SIMPLEX, dim, conn, coords);

  classify_elements(mesh, stk_mesh, stk_elems, merging_block_names);
  classify_sides(mesh, stk_mesh, sorted_nodes);
  finalize_classification(mesh);

  end_code();
}

void read_from_stk(Mesh* mesh,
                   const stk::mesh::BulkData & stk_mesh,
                   const std::vector<std::pair<std::string,std::set<std::string>>> & merging_block_names)
{
  std::vector<stk::mesh::Entity> stk_elems;
  stk::mesh::get_selected_entities( stk_mesh.mesh_meta_data().locally_owned_part(), stk_mesh.buckets(stk::topology::ELEMENT_RANK), stk_elems );
  read_from_stk(mesh, stk_mesh, stk_elems, merging_block_names);
}

}  // end namespace stk_mesh

}  // end namespace Omega_h
