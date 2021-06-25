#include <Omega_h_adapt.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_profile.hpp>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<double>("desired-num-elements");
  cmdline.add_arg<std::string>("output.vtk");
  auto const world = lib.world();
  if (!cmdline.parse_final(world, &argc, argv)) {
    return -1;
  }
  Omega_h::ScopedTimer scoped_timer("main");
  auto const world_size = world->size();
  auto const world_rank = world->rank();
  auto const inpath = cmdline.get<std::string>("input.osh");
  auto const desired_nelems = cmdline.get<double>("desired-num-elements");
  auto const outpath = cmdline.get<std::string>("output.vtk");
  auto const desired_nelems_per_rank = desired_nelems / world_size;
  Omega_h::Mesh mesh(&lib);
  for (int shift = 0; (1 << shift) <= world_size; ++shift) {
    int const group_size = (1 << shift);
    if (world_rank == 0)
      std::cout << "going to group size " << group_size << '\n';
    int const group_in_world = world_rank / group_size;
    int const group_rank = world_rank % group_size;
    auto const group = world->split(group_in_world, group_rank);
    OMEGA_H_CHECK(group->rank() == group_rank);
    if (group_in_world == 0) {
      if (shift == 0) {
        if (world_rank == 0) std::cout << "reading mesh...\n";
        mesh = Omega_h::binary::read(inpath, group, true);
      } else {
        if (world_rank == 0) std::cout << "repartitioning...\n";
        mesh.set_comm(group);
        mesh.balance();
      }
      Omega_h::AdaptOpts opts(&mesh);
      auto nelems = mesh.nglobal_ents(mesh.dim());
      if (world_rank == 0)
        std::cout << "mesh has " << nelems << " total elements\n";
      auto const desired_group_nelems = desired_nelems_per_rank * group_size;
      if (double(nelems) >= desired_group_nelems) {
        if (world_rank == 0)
          std::cout << "element count " << nelems << " >= target "
                    << desired_group_nelems << ", will not adapt\n";
      }

      auto vert_rc_ids = (mesh.ask_revClass(0)).ab2b;
      auto nbvert = vert_rc_ids.size();
      auto edge_rc_ids = (mesh.ask_revClass(1)).ab2b;
      auto nbedge = edge_rc_ids.size();
      auto face_rc_ids = (mesh.ask_revClass(2)).ab2b;
      auto nbface = face_rc_ids.size();
      auto reg_rc_ids = (mesh.ask_revClass(3)).ab2b;
      auto nbreg = reg_rc_ids.size();

      Omega_h::Write<Omega_h::Real> vals_v (nbvert, 100);
      Omega_h::Write<Omega_h::Real> vals_e (nbedge, 100);
      Omega_h::Write<Omega_h::Real> vals_f (nbface, 100);
      Omega_h::Write<Omega_h::Real> vals_r (nbreg, 100);
      mesh.add_rcField<Omega_h::Real> (0, "field", 1, Omega_h::Reals(vals_v));
      mesh.add_rcField<Omega_h::Real> (1, "field", 1, Omega_h::Reals(vals_e));
      mesh.add_rcField<Omega_h::Real> (2, "field", 1, Omega_h::Reals(vals_f));
      mesh.add_rcField<Omega_h::Real> (3, "field", 1, Omega_h::Reals(vals_r));
      add_rcField_transferMap (&opts, "field", OMEGA_H_INHERIT);

      while (double(nelems) < desired_group_nelems) {
        if (world_rank == 0)
          std::cout << "element count " << nelems << " < target "
                    << desired_group_nelems << ", will adapt\n";
        if (!mesh.has_tag(0, "metric")) {
          if (world_rank == 0)
            std::cout
                << "mesh had no metric, adding implied and adapting to it\n";
          Omega_h::add_implied_metric_tag(&mesh);
          Omega_h::adapt(&mesh, opts);
          nelems = mesh.nglobal_ents(mesh.dim());
          if (world_rank == 0)
            std::cout << "mesh now has " << nelems << " total elements\n";
        }
        auto metrics = mesh.get_array<double>(0, "metric");
        metrics = Omega_h::multiply_each_by(metrics, 1.2);
        auto const metric_ncomps =
            Omega_h::divide_no_remainder(metrics.size(), mesh.nverts());
        mesh.add_tag(0, "metric", metric_ncomps, metrics);
        if (world_rank == 0) std::cout << "adapting to scaled metric\n";
        Omega_h::adapt(&mesh, opts);
        nelems = mesh.nglobal_ents(mesh.dim());
        if (world_rank == 0)
          std::cout << "mesh now has " << nelems << " total elements\n";
      }
    }
    world->barrier();
  }
  Omega_h::vtk::FullWriter writer;
  writer = Omega_h::vtk::FullWriter(outpath, &mesh);
  writer.write();
  return 0;
}
