#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_array_ops.hpp>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  auto cmdline = Omega_h::CmdLine();
  cmdline.add_arg<std::string>("input.meshb");
  cmdline.add_arg<double>("target-num-elements");
  if (!cmdline.parse_final(world, &argc, argv)) {
    return -1;
  }
  auto target_nelems = cmdline.get<double>("target-num-elements");
  auto nelems_per_rank = target_nelems / world->size();
  auto world_rank = world->rank();
  auto self_comm = lib.self();
  Omega_h::Mesh mesh(&lib);
  if (world_rank == 0) {
    std::cout << "targeting " << target_nelems << " total elements, " << nelems_per_rank << " elements per rank\n";
    auto inmeshpath = cmdline.get<std::string>("input.meshb");
    Omega_h::meshb::read(&mesh, inmeshpath);
    std::cout << "rank 0 read " << inmeshpath << '\n';
    Omega_h::vtk::write_vtu("input.vtu", &mesh);
    std::cout << "Currently have exactly " << mesh.nelems() << " elements\n";
    std::cout << "getting implied metric...\n";
    auto implied_metric = Omega_h::get_implied_metrics(&mesh);
    auto predicted_nelems = Omega_h::get_expected_nelems(&mesh, implied_metric);
    std::cout << "Estimated " << predicted_nelems << " for the implied metric\n";
    mesh.add_tag(0, "metric", Omega_h::symm_ncomps(mesh.dim()), implied_metric);
    auto adapt_opts = Omega_h::AdaptOpts(&mesh);
    adapt_opts.verbosity = Omega_h::EXTRA_STATS;
    std::cout << "Adapting the mesh to fix quality\n";
    while (mesh.min_quality() < 0.2) {
      Omega_h::adapt(&mesh, adapt_opts);
    }
    implied_metric = mesh.get_array<Omega_h::Real>(0, "metric");
    predicted_nelems = Omega_h::get_expected_nelems(&mesh, implied_metric);
    std::cout << "Estimated " << predicted_nelems << " elements for the implied metric\n";
    std::cout << "Limiting gradation of the implied metric...\n";
    auto graded_implied_metric = Omega_h::limit_metric_gradation(&mesh, implied_metric, 1.0);
    std::cout << "Scaling the graded implied metric to get " << nelems_per_rank << " elements on rank 0\n";
    auto scalar = Omega_h::get_metric_scalar_for_nelems(&mesh, graded_implied_metric, nelems_per_rank);
    std::cout << "Metric scalar was " << scalar << '\n';
    auto target_metric = Omega_h::multiply_each_by(graded_implied_metric, scalar);
    predicted_nelems = Omega_h::get_expected_nelems(&mesh, target_metric);
    std::cout << "Estimated " << predicted_nelems << " elements for the target metric\n";
    mesh.add_tag(0, "target_metric", Omega_h::symm_ncomps(mesh.dim()), target_metric);
    std::cout << "Adapting the mesh to get " << nelems_per_rank << " elements on rank 0\n";
    while (Omega_h::approach_metric(&mesh, adapt_opts)) {
      Omega_h::adapt(&mesh, adapt_opts);
    }
    Omega_h::adapt(&mesh, adapt_opts);
    Omega_h::vtk::write_vtu("rank0-final.vtu", &mesh);
  }
  world->barrier();
  for (int subcomm_size = 2; subcomm_size <= world->size(); subcomm_size *= 2) {
    if (world_rank == 0) std::cout << "Moving to " << subcomm_size << " MPI ranks\n";
    auto is_in_subcomm = world_rank < subcomm_size;
    auto subcomm = world->split(int(is_in_subcomm), 0);
    if (is_in_subcomm) {
      mesh.set_comm(subcomm);
      subcomm->barrier();
      if (world_rank == 0) std::cout << "Balancing...\n";
      mesh.balance();
      subcomm->barrier();
      if (world_rank == 0) std::cout << "Ghosting...\n";
      mesh.set_parting(OMEGA_H_GHOSTED);
      subcomm->barrier();
      if (world_rank == 0) std::cout << "Current element count: " << mesh.nglobal_ents(mesh.dim()) << '\n';
      auto implied_metric = mesh.get_array<Omega_h::Real>(0, "metric");
      subcomm->barrier();
      auto predicted_nelems = Omega_h::get_expected_nelems(&mesh, implied_metric);
      subcomm->barrier();
      if (world_rank == 0) std::cout << "Estimated " << predicted_nelems << " elements for the implied metric\n";
      auto nsubcomm_elems = nelems_per_rank * subcomm_size;
      subcomm->barrier();
      if (world_rank == 0) std::cout << "Scaling the implied metric to get " << nsubcomm_elems << " elements on rank 0\n";
      auto scalar = Omega_h::get_metric_scalar_for_nelems(&mesh, implied_metric, nsubcomm_elems);
      if (world_rank == 0) std::cout << "Metric scalar was " << scalar << '\n';
      auto target_metric = Omega_h::multiply_each_by(implied_metric, scalar);
      predicted_nelems = Omega_h::get_expected_nelems(&mesh, target_metric);
      if (world_rank == 0) std::cout << "Estimated " << predicted_nelems << " for the target metric\n";
      mesh.add_tag(0, "target_metric", Omega_h::symm_ncomps(mesh.dim()), target_metric);
      if (world_rank == 0) std::cout << "Adapting the mesh to get " << nelems_per_rank << " elements on rank 0\n";
      auto adapt_opts = Omega_h::AdaptOpts(&mesh);
      adapt_opts.verbosity = Omega_h::EXTRA_STATS;
      while (Omega_h::approach_metric(&mesh, adapt_opts)) {
        Omega_h::adapt(&mesh, adapt_opts);
      }
      Omega_h::adapt(&mesh, adapt_opts);
    }
    world->barrier();
  }
}
