#include <Omega_h_file.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_profile.hpp>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<double>("desired-num-elements");
  cmdline.add_arg<std::string>("output.osh");
  auto const world = lib.world();
  if (!cmdline.parse_final(world, &argc, argv)) {
    return -1;
  }
  Omega_h::ScopedTimer scoped_timer("main");
  auto const world_size = world->size();
  auto const world_rank = world->rank();
  auto const inpath = cmdline.get<std::string>("input.osh");
  auto const outpath = cmdline.get<std::string>("output.osh");
  Omega_h::Mesh mesh(&lib);
  for (int shift = 0; (1 << shift) <= world_size; ++shift) {
    int const group_size = (1 << shift);
    int const group_in_world = world_rank / group_size;
    int const group_rank = world_rank % group_size;
    auto const group = world->split(group_in_world, group_rank);
    OMEGA_H_CHECK(group->rank() == group_rank);
    if (group_in_world == 0) {
      if (shift == 0) {
        std::cout << "reading mesh...\n";
        mesh = Omega_h::binary::read(inpath, group, true);
      } else {
        mesh.set_comm(group);
        mesh.balance();
      }
      auto const nelems = mesh.nglobal_ents(mesh.dim());
      auto const nelems_local = mesh.nelems();
      if (world_rank == 0) std::cout << "mesh has " << nelems << " total elements, " << nelems_local << " on rank 0\n";
    }
    world->barrier();
  }
  return 0;
}
