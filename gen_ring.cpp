#include <gmodel.hpp>

int main() {
  gmod::default_size = 0.5;
  auto disk = new_disk(gmod::Vector{0, 0, 0}, gmod::Vector{0, 0, 1},
                       gmod::Vector{4, 0, 0});
  auto inner = new_circle(gmod::Vector{0, 0, 0}, gmod::Vector{0, 0, 1},
                          gmod::Vector{2, 0, 0});
  gmod::add_use(disk, gmod::REVERSE, inner);
  gmod::write_closure_to_geo(disk, "ring.geo");
}
