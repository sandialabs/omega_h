#include <gmodel.hpp>

int main() {
  auto b = gmod::new_ball(gmod::Vector{0.5, 0.5, 0.5}, gmod::Vector{0, 0, 1},
                          gmod::Vector{0.2, 0, 0});
  gmod::write_closure_to_geo(b, "ball.geo");
}
