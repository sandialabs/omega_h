#include <gmodel.hpp>

int main() {
  auto g = gmod::new_group();
  gmod::default_size = 0.1;
  auto fluid = new_cube(gmod::Vector{0, 0, 0}, gmod::Vector{1, 0, 0},
                        gmod::Vector{0, 1, 0}, gmod::Vector{0, 0, 1.5});
  gmod::add_to_group(g, fluid);
  gmod::default_size = 0.02;
  auto cyl_disk = new_disk(gmod::Vector{0.5, 0.5, 0.1}, gmod::Vector{0, 0, 1},
                           gmod::Vector{1.0 / 8.0, 0, 0});
  auto cylinder = gmod::extrude_face(cyl_disk, gmod::Vector{0, 0, 0.5}).middle;
  gmod::add_to_group(g, cylinder);
  gmod::insert_into(fluid, cylinder);
  auto tube_disk = new_disk(gmod::Vector{0.5, 0.5, 0.7}, gmod::Vector{0, 0, 1},
                            gmod::Vector{1.0 / 8.0 + 2.0 / 32.0, 0, 0});
  auto tube_inner =
      new_circle(gmod::Vector{0.5, 0.5, 0.7}, gmod::Vector{0, 0, 1},
                 gmod::Vector{1.0 / 8.0 + 1.0 / 32.0, 0, 0});
  gmod::add_use(tube_disk, gmod::REVERSE, tube_inner);
  auto tube = extrude_face(tube_disk, gmod::Vector{0, 0, 0.5}).middle;
  gmod::add_to_group(g, tube);
  gmod::insert_into(fluid, tube);
  gmod::write_closure_to_geo(g, "cylinder_thru_tube.geo");
}
