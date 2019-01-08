#include <Omega_h_compare.hpp>
#include <Omega_h_vtk.hpp>

#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  CmdLine cmdline;
  get_diff_program_cmdline("gold_dir", "result_dir", &cmdline);
  cmdline.add_flag("-onestep", "Each directory is one time step");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) {
    return -1;
  }
  auto patha = cmdline.get<std::string>("gold_dir");
  auto pathb = cmdline.get<std::string>("result_dir");
  MeshCompareOpts opts;
  Omega_h_Comparison max_result = OMEGA_H_SAME;
  auto one_step = cmdline.parsed("-onestep");
  if (one_step) {
    auto pvtupatha = vtk::get_pvtu_path(patha);
    auto pvtupathb = vtk::get_pvtu_path(pathb);
    Mesh mesha(&lib);
    vtk::read_parallel(pvtupatha, lib.world(), &mesha);
    accept_diff_program_cmdline(cmdline, &mesha, &opts, &max_result);
    Mesh meshb(&lib);
    vtk::read_parallel(pvtupathb, lib.world(), &meshb);
    auto res = compare_meshes(&mesha, &meshb, opts, true, false);
    if (res <= max_result) return 0;
    return 2;
  }
  auto const pvdpatha = vtk::get_pvd_path(patha);
  auto const pvdpathb = vtk::get_pvd_path(pathb);
  std::vector<Real> timesa;
  std::vector<filesystem::path> pvtupathsa;
  vtk::read_pvd(pvdpatha, &timesa, &pvtupathsa);
  std::vector<Real> timesb;
  std::vector<filesystem::path> pvtupathsb;
  vtk::read_pvd(pvdpathb, &timesb, &pvtupathsb);
  if (timesa.size() != timesb.size()) {
    std::cout << "different number of time steps (" << timesa.size() << " vs "
              << timesb.size() << ")\n";
    return 2;
  }
  for (std::size_t step = 0; step < timesa.size(); ++step) {
    auto pvtupatha = pvtupathsa[step];
    auto pvtupathb = pvtupathsb[step];
    Mesh mesha(&lib);
    vtk::read_parallel(pvtupatha, lib.world(), &mesha);
    if (step == 0)
      accept_diff_program_cmdline(cmdline, &mesha, &opts, &max_result);
    auto timea = timesa[step];
    auto timeb = timesb[step];
    if (!compare_real(timea, timeb, opts.time_step_opts)) {
      std::cout << "time for step " << step << " differs (" << timea << " vs "
                << timeb << ")\n";
      return 2;
    }
    Mesh meshb(&lib);
    vtk::read_parallel(pvtupathb, lib.world(), &meshb);
    auto res = compare_meshes(&mesha, &meshb, opts, true, false);
    if (res == OMEGA_H_DIFF) {
      std::cout << "step " << step << " differs.\n";
      return 2;
    }
    if (res == OMEGA_H_MORE && max_result < OMEGA_H_MORE) {
      std::cout << "step " << step
                << " is a superset but \"-superset\" not given\n";
      return 2;
    }
  }
  return 0;
}
