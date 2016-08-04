#include <gmodel.hpp>
#include <iostream>

int main() {
  double major = 0.15;
  double minor = 0.08;
  double blade_length = 0.3;
  double link_length = 0.2;
  double diameter = 0.05;
  double width = 1.7;
  double height = 1.2;
  double thickness = 0.4;
  double spacing = width - height;
  gmod::default_size = 0.005;
  auto edisk0 = gmod::new_elliptical_disk(gmod::Vector{0, 0, 0},
      gmod::Vector{major / 2, 0, 0},
      gmod::Vector{0, 0, minor / 2});
  auto edisk1 = gmod::copy_closure(edisk0);
  auto disk0 = gmod::new_disk(gmod::Vector{0, 0, 0},
      gmod::Vector{0, 1, 0}, gmod::Vector{diameter / 2, 0, 0});
  auto blade0 =
    gmod::extrude_face(edisk0, gmod::Vector{0, blade_length, 0}).middle;
  auto blade1 =
    gmod::extrude_face(edisk1, gmod::Vector{0, blade_length, 0}).middle;
  gmod::transform_closure(blade0, gmod::identity_matrix(),
      gmod::Vector{0, link_length / 2, 0});
  gmod::transform_closure(blade1,
      gmod::rotation_matrix(gmod::Vector{0, 0, 1}, gmod::PI),
      gmod::Vector{0, -link_length / 2, 0});
  auto link_extrusion = gmod::extrude_face(disk0, gmod::Vector{0, link_length, 0});
  auto link = link_extrusion.middle;
  auto disk1 = link_extrusion.end;
  gmod::transform_closure(link, gmod::identity_matrix(),
      gmod::Vector{0, -link_length / 2, 0});
  gmod::weld_volume_face_into(blade0, link, edisk0, disk1);
  gmod::weld_volume_face_into(blade1, link, edisk1, disk0);
  auto assembly0 = gmod::new_group();
  gmod::add_to_group(assembly0, blade0);
  gmod::add_to_group(assembly0, blade1);
  gmod::add_to_group(assembly0, link);
  auto assembly1 = gmod::copy_closure(assembly0);
  gmod::transform_closure(assembly0, gmod::identity_matrix(),
      gmod::Vector{-spacing / 2, 0, 0});
  gmod::transform_closure(assembly1,
      gmod::rotation_matrix(gmod::Vector{0, 0, 1}, gmod::PI / 2),
      gmod::Vector{spacing / 2, 0, 0});
  gmod::default_size = 0.04;
  auto fluid = gmod::new_cube(
      gmod::Vector{-width / 2, -height / 2, - thickness / 2},
      gmod::Vector{width, 0, 0},
      gmod::Vector{0, height, 0},
      gmod::Vector{0, 0, thickness});
  gmod::insert_into(fluid, assembly0);
  gmod::insert_into(fluid, assembly1);
  auto model = gmod::new_group();
  gmod::add_to_group(model, assembly0);
  gmod::add_to_group(model, assembly1);
  gmod::add_to_group(model, fluid);
  gmod::write_closure_to_geo(model, "twin_rotor.geo");
  std::cout << "assembly0 regions:";
  for (auto vol_use : assembly0->used)
    std::cout << ' ' << vol_use.obj->id;
  std::cout << '\n';
  std::cout << "assembly1 regions:";
  for (auto vol_use : assembly1->used)
    std::cout << ' ' << vol_use.obj->id;
  std::cout << '\n';
}
