#include <gmodel.hpp>

int main() {
  using V = gmod::Vector;
  gmod::default_size = 1./18.;
  auto g = gmod::new_group();
  auto air = gmod::new_square(V{0, 0, 0}, V{1, 0, 0}, V{0, 1, 0});
  auto weight = gmod::new_square(
      V{1./6., 3./8., 0}, V{1./4., 0, 0}, V{0, 1./4., 0});
  auto droplet = gmod::new_disk(
      V{17./24., 1./2., 0}, V{0, 0, 1}, V{1./8., 0, 0});
  gmod::insert_into(air, weight);
  gmod::insert_into(air, droplet);
  gmod::add_to_group(g, air);
  gmod::add_to_group(g, weight);
  gmod::add_to_group(g, droplet);
  gmod::write_closure_to_geo(g, "square_and_disk_in_square.geo");
  printf("weight is %d\n", weight->id);
  printf("droplet is %d\n", droplet->id);
}
