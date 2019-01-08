#include "Omega_h_compare.hpp"

#include <iostream>

#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

bool check_regression(filesystem::path const& prefix, Mesh* mesh) {
  auto comm = mesh->comm();
  auto goldpath = prefix;
  goldpath += ".osh";
  if (!exists(goldpath)) {
    if (comm->size() == 1) {
      binary::write(goldpath, mesh);
      if (!mesh->library()->silent_) {
        std::cout << "gold path \"" << goldpath << "\" did not exist yet,\n";
        std::cout << "created it from this run.\n";
      }
      return true;
    } else {
      auto tmppath = prefix;
      tmppath += "_tmp.osh";
      binary::write(tmppath, mesh);
      if (!mesh->library()->silent_ && comm->rank() == 0) {
        std::cout << "gold path \"" << goldpath;
        std::cout << "\" does not exist and this run is parallel.\n";
        std::cout << "If you really want to use this run as the gold, do:\n";
        std::cout << "  mv \"" << tmppath << "\" \"" << goldpath << "\"\n";
      }
      return false;
    }
  }
  Mesh gold_mesh(mesh->library());
  binary::read(goldpath, comm, &gold_mesh);
  auto opts = MeshCompareOpts::init(mesh, VarCompareOpts::zero_tolerance());
  auto res = compare_meshes(&gold_mesh, mesh, opts, true);
  if (res == OMEGA_H_SAME) {
    if (!mesh->library()->silent_ && comm->rank() == 0) {
      std::cout << "This run matches gold \"" << goldpath << "\"\n";
    }
    return true;
  }
  if (res == OMEGA_H_MORE) {
    auto newpath = prefix;
    newpath += "_new.osh";
    binary::write(newpath, mesh);
    if (!mesh->library()->silent_ && comm->rank() == 0) {
      std::cout << "This run, stored at \"" << newpath << "\",\n";
      std::cout << "has more tags than \"" << goldpath << "\"\n";
      std::cout << "It should probably be made the new gold, like this:\n";
      std::cout << "  rm -rf \"" << goldpath << "\"\n";
      std::cout << "  mv \"" << newpath << "\" \"" << goldpath << "\"\n";
    }
    return false;
  }
  auto badpath = prefix;
  badpath += "_bad.osh";
  binary::write(badpath, mesh);
  if (!mesh->library()->silent_ && comm->rank() == 0) {
    std::cout << "This run, stored at \"" << badpath << "\",\n";
    std::cout << "does not match the gold at \"" << goldpath << "\"\n";
  }
  return false;
}

}  // end namespace Omega_h
