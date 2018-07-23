#include <Omega_h_adapt.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto t0 = Omega_h::now();
  Omega_h::Now t1;
  auto world = lib.world();
  auto cmdline = Omega_h::CmdLine();
  cmdline.add_arg<std::string>("input.meshb");
  cmdline.add_arg<double>("final-complexity");
  if (!cmdline.parse_final(world, &argc, argv)) {
    return -1;
  }
  auto final_complexity = cmdline.get<double>("final-complexity");
  auto final_nelems =
      Omega_h::get_expected_nelems_from_complexity(final_complexity, 3);
  auto nelems_per_rank = final_nelems / world->size();
  auto world_rank = world->rank();
  if (world_rank == 0)
    std::cout << "Estimating " << final_nelems
              << " final elements for the final complexity of "
              << final_complexity << '\n';
  auto self_comm = lib.self();
  Omega_h::Mesh mesh(&lib);
  if (world_rank == 0) {
    std::cout << "targeting " << final_nelems << " total elements, "
              << nelems_per_rank << " elements per rank\n";
    auto inmeshpath = cmdline.get<std::string>("input.meshb");
    Omega_h::meshb::read(&mesh, inmeshpath);
    std::cout << "rank 0 read " << inmeshpath << '\n';
    // Omega_h::vtk::write_vtu("input.vtu", &mesh);
    std::cout << "Currently have exactly " << mesh.nelems() << " elements\n";
    std::cout << "getting implied metric...\n";
    auto implied_metric = Omega_h::get_implied_metrics(&mesh);
    auto complexity = Omega_h::get_complexity(&mesh, implied_metric);
    std::cout << "Initial complexity was " << complexity << '\n';
    auto predicted_nelems =
        Omega_h::get_expected_nelems_from_complexity(complexity, mesh.dim());
    std::cout << "Estimated " << predicted_nelems
              << " for the implied metric\n";
    mesh.add_tag(0, "metric", Omega_h::symm_ncomps(mesh.dim()), implied_metric);
    auto adapt_opts = Omega_h::AdaptOpts(&mesh);
    adapt_opts.verbosity = Omega_h::EXTRA_STATS;
    std::cout << "Adapting the mesh to fix quality\n";
    while (mesh.min_quality() < 0.2) {
      Omega_h::adapt(&mesh, adapt_opts);
    }
    t1 = Omega_h::now();
    implied_metric = mesh.get_array<Omega_h::Real>(0, "metric");
    std::cout << "Limiting gradation of the implied metric...\n";
    implied_metric =
        Omega_h::limit_metric_gradation(&mesh, implied_metric, 1.0);
    std::cout << "Scaling the graded implied metric to get " << nelems_per_rank
              << " elements on rank 0\n";
    auto scalar = Omega_h::get_metric_scalar_for_nelems(
        mesh.dim(), mesh.nelems(), nelems_per_rank);
    std::cout << "Metric scalar was " << scalar << '\n';
    auto target_metric = Omega_h::multiply_each_by(implied_metric, scalar);
    mesh.add_tag(
        0, "target_metric", Omega_h::symm_ncomps(mesh.dim()), target_metric);
    std::cout << "Adapting the mesh to satisfy the target metric...\n";
    while (Omega_h::approach_metric(&mesh, adapt_opts)) {
      Omega_h::adapt(&mesh, adapt_opts);
    }
    Omega_h::adapt(&mesh, adapt_opts);
    std::cout << "wound up with " << mesh.nelems() << " elements\n";
    // Omega_h::vtk::write_vtu("rank0-final.vtu", &mesh);
  }
  world->barrier();
  for (int subcomm_size = 2; subcomm_size <= world->size(); subcomm_size *= 2) {
    if (world_rank == 0)
      std::cout << "Moving to " << subcomm_size << " MPI ranks\n";
    auto is_in_subcomm = world_rank < subcomm_size;
    auto subcomm = world->split(int(is_in_subcomm), 0);
    if (is_in_subcomm) {
      mesh.set_comm(subcomm);
      if (world_rank == 0) std::cout << "Balancing...\n";
      mesh.balance();
      // Omega_h::vtk::write_parallel("balanced", &mesh,
      // Omega_h::vtk::dont_compress);
      if (world_rank == 0) std::cout << "Ghosting...\n";
      mesh.set_parting(OMEGA_H_GHOSTED);
      // Omega_h::vtk::write_parallel("ghosted", &mesh,
      // Omega_h::vtk::dont_compress);
      auto nglobal_elems = mesh.nglobal_ents(mesh.dim());
      if (world_rank == 0)
        std::cout << "Current element count: " << nglobal_elems << '\n';
      auto implied_metric = mesh.get_array<Omega_h::Real>(0, "metric");
      auto nsubcomm_elems = nelems_per_rank * subcomm_size;
      if (world_rank == 0)
        std::cout << "Scaling the implied metric to get " << nsubcomm_elems
                  << " elements in " << subcomm_size << " ranks\n";
      auto scalar = Omega_h::get_metric_scalar_for_nelems(
          mesh.dim(), nglobal_elems, nsubcomm_elems);
      if (world_rank == 0) std::cout << "Metric scalar was " << scalar << '\n';
      auto target_metric = Omega_h::multiply_each_by(implied_metric, scalar);
      mesh.add_tag(
          0, "target_metric", Omega_h::symm_ncomps(mesh.dim()), target_metric);
      if (world_rank == 0)
        std::cout << "Adapting the mesh to get " << nsubcomm_elems
                  << " elements in " << subcomm_size << " ranks\n";
      auto adapt_opts = Omega_h::AdaptOpts(&mesh);
      adapt_opts.verbosity = Omega_h::EXTRA_STATS;
      while (Omega_h::approach_metric(&mesh, adapt_opts)) {
        Omega_h::adapt(&mesh, adapt_opts);
      }
      Omega_h::adapt(&mesh, adapt_opts);
      nglobal_elems = mesh.nglobal_ents(mesh.dim());
      if (world_rank == 0)
        std::cout << "Wound up with " << nglobal_elems << " elements\n";
    }
    world->barrier();
  }
  auto t2 = Omega_h::now();
  if (world_rank == 0)
    std::cout << "Tprep " << (t1 - t0) << " Tadapt " << (t2 - t1) << " Ttotal "
              << (t2 - t0) << '\n';
}
