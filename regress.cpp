bool check_regression(std::string const& prefix, Mesh* mesh,
    Real tol, Real floor) {
  auto comm = mesh->comm();
  auto goldpath = prefix + ".osh";
  if (!directory_exists(goldpath.c_str())) {
    if (comm->size() == 1) {
      binary::write(goldpath, mesh);
      std::cout << "gold path \"" << goldpath << "\" did not exist yet,\n";
      std::cout << "created it from this run.\n";
      return true;
    } else if (comm->rank() == 0) {
      auto tmppath = prefix + "_tmp.osh";
      binary::write(tmppath, mesh);
      std::cout << "gold path \"" << goldpath;
      std::cout << "\" does not exist and this run is parallel.\n";
      std::cout << "If you really want to use this run as the gold, do:\n";
      std::cout << "  mv \"" << tmppath << "\" \"" << goldpath << "\"\n";
      return false;
    }
  }
  Mesh mesh2;
  binary::read(goldpath, comm, &mesh2);
  if (compare_meshes(mesh, &mesh2, tol, floor, false, true)) {
    if (comm->rank() == 0) {
      std::cout << "This run matches gold \"" << goldpath << "\"\n";
    }
    return true;
  }
  std::cout << "meshes not exactly the same, trying superset...\n";
  if (compare_meshes(mesh, &mesh2, tol, floor, true, true)) {
    auto newpath = prefix + "_new.osh";
    binary::write(newpath, mesh);
    if (comm->rank() == 0) {
      std::cout << "This run, stored at \"" << newpath << "\",\n";
      std::cout << "has more tags than \"" << goldpath << "\"\n";
      std::cout << "It should probably be made the new gold, like this:\n";
      std::cout << "  rm -rf \"" << goldpath << "\"\n",
      std::cout << "  mv \"" << newpath << "\" \"" << goldpath << "\"\n";
    }
    return true;
  }
  auto badpath = prefix + "_bad.osh";
  binary::write(badpath, mesh);
  if (comm->rank() == 0) {
    std::cout << "This run, stored at \"" << badpath << "\",\n";
    std::cout << "does not match the gold at \"" << goldpath << "\"\n";
  }
  return false;
}
