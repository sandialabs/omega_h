#include "Omega_h.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "algebra.hpp"
#include "vtk.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  Real tol = 1e-6;
  Real floor = 0.0;
  bool get_tol = false;
  bool get_floor = false;
  bool get_help = false;
  char const* patha = nullptr;
  char const* pathb = nullptr;
  bool allow_superset = false;
  bool one_step = false;
  for (int i = 1; i < argc; ++i) {
    if (get_tol) {
      tol = atof(argv[i]);
      get_tol = false;
      continue;
    }
    if (get_floor) {
      floor = atof(argv[i]);
      get_floor = false;
      continue;
    }
    if (!strcmp(argv[i], "-tolerance")) {
      get_tol = true;
      continue;
    }
    if (!strcmp(argv[i], "-Floor")) {
      get_floor = true;
      continue;
    }
    if (!strcmp(argv[i], "-help")) {
      get_help = true;
      continue;
    }
    if (!strcmp(argv[i], "-superset")) {
      allow_superset = true;
      continue;
    }
    if (!strcmp(argv[i], "-onestep")) {
      one_step = true;
      continue;
    }
    if (!patha) {
      patha = argv[i];
      continue;
    }
    if (!pathb) {
      pathb = argv[i];
      continue;
    }
  }
  if (get_tol) {
    printf("-tolerance needs an argument\n");
    return -1;
  }
  if (get_floor) {
    printf("-Floor needs an argument\n");
    return -1;
  }
  if (get_help || !patha || !pathb) {
    std::cout << "\n";
    std::cout << "usage: vtkdiff [options] gold_dir result_dir\n";
    std::cout << "   or: vtkdiff [-help]               (usage)\n";
    std::cout << "\n";
    std::cout << "    -help (Print this summary and exit.)\n";
    std::cout << "    -tolerance <$val> (Overrides the default tolerance of "
                 "1.0E-6.)\n";
    std::cout << "    -Floor <$val> (Overrides the default floor tolerance of "
                 "0.0.)\n";
    std::cout << "    -superset (Allow result to have more arrays than gold)\n";
    std::cout << "    -onestep (Expect only a .pvtu for a single time step)\n";
    return -1;
  }
  if (one_step) {
    auto pvtupatha = vtk::get_pvtu_path(patha);
    auto pvtupathb = vtk::get_pvtu_path(pathb);
    Mesh mesha;
    vtk::read_parallel(pvtupatha, lib.world(), &mesha);
    Mesh meshb;
    vtk::read_parallel(pvtupathb, lib.world(), &meshb);
    auto res = compare_meshes(&mesha, &meshb, tol, floor, true, false);
    if (res == OMEGA_H_SAME) return 0;
    if (allow_superset && res == OMEGA_H_MORE) return 0;
    return 2;
  }
  auto pvdpatha = vtk::get_pvd_path(patha);
  auto pvdpathb = vtk::get_pvd_path(pathb);
  std::vector<Real> timesa;
  std::vector<std::string> pvtupathsa;
  vtk::read_pvd(pvdpatha, &timesa, &pvtupathsa);
  std::vector<Real> timesb;
  std::vector<std::string> pvtupathsb;
  vtk::read_pvd(pvdpathb, &timesb, &pvtupathsb);
  if (timesa.size() != timesb.size()) {
    std::cout << "different number of time steps (" << timesa.size() << " vs "
              << timesb.size() << ")\n";
    return 2;
  }
  for (std::size_t step = 0; step < timesa.size(); ++step) {
    auto timea = timesa[step];
    auto timeb = timesb[step];
    if (!are_close(timea, timeb, tol, floor)) {
      std::cout << "time for step " << step << " differs (" << timea << " vs "
                << timeb << ")\n";
      return 2;
    }
    auto pvtupatha = pvtupathsa[step];
    auto pvtupathb = pvtupathsb[step];
    Mesh mesha;
    vtk::read_parallel(pvtupatha, lib.world(), &mesha);
    Mesh meshb;
    vtk::read_parallel(pvtupathb, lib.world(), &meshb);
    auto res = compare_meshes(&mesha, &meshb, tol, floor, true, false);
    if (res == OMEGA_H_DIFF) {
      std::cout << "step " << step << " differs.\n";
      return 2;
    }
    if (res == OMEGA_H_MORE && !allow_superset) {
      std::cout << "step " << step
                << " is a superset but \"-superset\" not given\n";
      return 2;
    }
  }
  return 0;
}
