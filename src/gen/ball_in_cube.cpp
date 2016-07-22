#include <gmodel.hpp>

int main() {
  auto g = gmod::new_group();
  auto fluid = gmod::new_cube(gmod::Vector{0, 0, 0}, gmod::Vector{1, 0, 0},
                              gmod::Vector{0, 1, 0}, gmod::Vector{0, 0, 1});
  gmod::add_to_group(g, fluid);
  auto solid = gmod::new_ball(gmod::Vector{0.5, 0.5, 0.5},
                              gmod::Vector{0, 0, 1}, gmod::Vector{0.2, 0, 0});
  gmod::add_to_group(g, solid);
  gmod::insert_into(fluid, solid);
  gmod::write_closure_to_geo(g, "ball_in_cube.geo");
}
