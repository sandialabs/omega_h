#include <gmodel.hpp>

int main() {
  auto g = gmod::new_group();
  auto fluid = gmod::new_cube(gmod::Vector{0, 0, 0}, gmod::Vector{1, 0, 0},
                              gmod::Vector{0, 1, 0}, gmod::Vector{0, 0, 2});
  gmod::add_to_group(g, fluid);
  auto solid1 = gmod::new_ball(gmod::Vector{0.5, 0.5, 0.5},
                               gmod::Vector{0, 0, 1}, gmod::Vector{0.2, 0, 0});
  gmod::add_to_group(g, solid1);
  gmod::insert_into(fluid, solid1);
  auto solid2 = gmod::new_ball(gmod::Vector{0.5, 0.5, 1.5},
                               gmod::Vector{0, 0, 1}, gmod::Vector{0.2, 0, 0});
  gmod::add_to_group(g, solid2);
  gmod::insert_into(fluid, solid2);
  gmod::write_closure_to_geo(g, "balls_in_box.geo");
  printf("solids are %u and %u\n", solid1->id, solid2->id);
}
