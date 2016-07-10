#include "internal.hpp"

#include "file.hpp"

namespace osh {

bool check_regression(std::string const& prefix, Mesh* mesh, Real tol,
                      Real floor) {
  auto comm = mesh->comm();
  auto goldpath = prefix + ".osh";
  if (!directory_exists(goldpath.c_str())) {
    if (comm->size() == 1) {
      binary::write(goldpath, mesh);
      std::cout << "gold path \"" << goldpath << "\" did not exist yet,\n";
      std::cout << "created it from this run.\n";
      return true;
    } else {
      auto tmppath = prefix + "_tmp.osh";
      binary::write(tmppath, mesh);
      if (comm->rank() == 0) {
        std::cout << "gold path \"" << goldpath;
        std::cout << "\" does not exist and this run is parallel.\n";
        std::cout << "If you really want to use this run as the gold, do:\n";
        std::cout << "  mv \"" << tmppath << "\" \"" << goldpath << "\"\n";
      }
      return false;
    }
  }
  Mesh gold_mesh;
  binary::read(goldpath, comm, &gold_mesh);
  auto res = compare_meshes(&gold_mesh, mesh, tol, floor, true);
  if (res == OSH_SAME) {
    if (comm->rank() == 0) {
      std::cout << "This run matches gold \"" << goldpath << "\"\n";
    }
    return true;
  }
  if (res == OSH_MORE) {
    auto newpath = prefix + "_new.osh";
    binary::write(newpath, mesh);
    if (comm->rank() == 0) {
      std::cout << "This run, stored at \"" << newpath << "\",\n";
      std::cout << "has more tags than \"" << goldpath << "\"\n";
      std::cout << "It should probably be made the new gold, like this:\n";
      std::cout << "  rm -rf \"" << goldpath << "\"\n";
      std::cout << "  mv \"" << newpath << "\" \"" << goldpath << "\"\n";
    }
    return false;
  }
  auto badpath = prefix + "_bad.osh";
  binary::write(badpath, mesh);
  if (comm->rank() == 0) {
    std::cout << "This run, stored at \"" << badpath << "\",\n";
    std::cout << "does not match the gold at \"" << goldpath << "\"\n";
  }
  return false;
}

} //end namespace osh
