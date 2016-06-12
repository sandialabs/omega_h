bool warp_to_limit(Mesh* mesh, Real min_qual) {
  if (!mesh.has_tag(VERT, "warp")) return false;
  CHECK(mesh.min_quality() >= min_qual);
  auto coords = mesh.coords();
  auto warp = mesh.get_array<Real>(VERT, "warp");
  mesh.set_coords(add_each(coords, warp));
  if (mesh.min_quality() >= min_qual) {
    mesh.remove_tag(VERT, "warp");
    return true;
  }
  auto remainder = Reals(warp.size(), 0.0);
  do {
    auto half_warp = multiply_each_by(1.0 / 2.0, warp);
    warp = half_warp;
    remainder = add_each(remainder, half_warp);
    mesh.set_coords(add_each(coords, warp));
  } while (mesh.min_quality() < min_qual);
  mesh.set_tag(VERT, "warp", remainder);
  return true;
}
