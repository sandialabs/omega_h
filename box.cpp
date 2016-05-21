void make_2d_box(Real x, Real y, LO nx, LO ny,
    LOs& qv2v, Reals& coords) {
  LO nq = nx * ny;
  LO nvx = nx + 1;
  LO nvy = ny + 1;
  LO nv = nvx * nvy;
  Real dx = x / nx;
  Real dy = y / ny;
  Write<Real> coords_(nv * 2);
  auto fill_coords = LAMBDA(LO v) {
    LO i = v % nvx;
    LO j = v / nvx;
    coords_[v * 2 + 0] = i * dx;
    coords_[v * 2 + 1] = j * dy;
  };
  parallel_for(nv, fill_coords);
  Write<LO> qv2v_(nq * 4);
  auto fill_conn = LAMBDA(LO q) {
    LO i = q % nx;
    LO j = q / nx;
    qv2v_[q * 4 + 0] = (j + 0) * nvx + (i + 0);
    qv2v_[q * 4 + 1] = (j + 0) * nvx + (i + 1);
    qv2v_[q * 4 + 2] = (j + 1) * nvx + (i + 1);
    qv2v_[q * 4 + 3] = (j + 1) * nvx + (i + 0);
  };
  parallel_for(nq, fill_conn);
  qv2v = qv2v_;
  coords = coords_;
}

void make_3d_box(
    Real x, Real y, Real z,
    LO nx, LO ny, LO nz,
    LOs& hv2v, Reals& coords) {
  LO nxy = nx * ny;
  LO nh = nx * ny * nz;
  LO nvx = nx + 1;
  LO nvy = ny + 1;
  LO nvz = nz + 1;
  LO nvxy = nvx * nvy;
  LO nv = nvx * nvy * nvz;
  Real dx = x / nx;
  Real dy = y / ny;
  Real dz = z / nz;
  Write<Real> coords_(nv * 3);
  auto fill_coords = LAMBDA(LO v) {
    LO ij = v % nvxy;
    LO k = v / nvxy;
    LO i = ij % nvx;
    LO j = ij / nvx;
    coords_[v * 3 + 0] = i * dx;
    coords_[v * 3 + 1] = j * dy;
    coords_[v * 3 + 2] = k * dz;
  };
  parallel_for(nv, fill_coords);
  Write<LO> hv2v_(nh * 8);
  auto fill_conn = LAMBDA(LO h) {
    LO ij = h % nxy;
    LO k = h / nxy;
    LO i = ij % nx;
    LO j = ij / nx;
    hv2v_[h * 8 + 0] = (k + 0) * nvxy + (j + 0) * nvx + (i + 0);
    hv2v_[h * 8 + 1] = (k + 0) * nvxy + (j + 0) * nvx + (i + 1);
    hv2v_[h * 8 + 2] = (k + 0) * nvxy + (j + 1) * nvx + (i + 1);
    hv2v_[h * 8 + 3] = (k + 0) * nvxy + (j + 1) * nvx + (i + 0);
    hv2v_[h * 8 + 4] = (k + 1) * nvxy + (j + 0) * nvx + (i + 0);
    hv2v_[h * 8 + 5] = (k + 1) * nvxy + (j + 0) * nvx + (i + 1);
    hv2v_[h * 8 + 6] = (k + 1) * nvxy + (j + 1) * nvx + (i + 1);
    hv2v_[h * 8 + 7] = (k + 1) * nvxy + (j + 1) * nvx + (i + 0);
  };
  parallel_for(nh, fill_conn);
  hv2v = hv2v_;
  coords = coords_;
}
