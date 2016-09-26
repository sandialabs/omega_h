#include "box.hpp"

#include "algebra.hpp"
#include "access.hpp"
#include "loop.hpp"

namespace Omega_h {

void make_2d_box(
    Real x, Real y, LO nx, LO ny, LOs* qv2v_out, Reals* coords_out) {
  LO nq = nx * ny;
  LO nvx = nx + 1;
  LO nvy = ny + 1;
  LO nv = nvx * nvy;
  Real dx = x / nx;
  Real dy = y / ny;
  Write<Real> coords(nv * 2);
  auto fill_coords = LAMBDA(LO v) {
    LO i = v % nvx;
    LO j = v / nvx;
    coords[v * 2 + 0] = i * dx;
    coords[v * 2 + 1] = j * dy;
  };
  parallel_for(nv, fill_coords);
  Write<LO> qv2v(nq * 4);
  auto fill_conn = LAMBDA(LO q) {
    LO i = q % nx;
    LO j = q / nx;
    qv2v[q * 4 + 0] = (j + 0) * nvx + (i + 0);
    qv2v[q * 4 + 1] = (j + 0) * nvx + (i + 1);
    qv2v[q * 4 + 2] = (j + 1) * nvx + (i + 1);
    qv2v[q * 4 + 3] = (j + 1) * nvx + (i + 0);
  };
  parallel_for(nq, fill_conn);
  *qv2v_out = qv2v;
  *coords_out = coords;
}

void make_3d_box(Real x, Real y, Real z, LO nx, LO ny, LO nz, LOs* hv2v_out,
    Reals* coords_out) {
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
  Write<Real> coords(nv * 3);
  auto fill_coords = LAMBDA(LO v) {
    LO ij = v % nvxy;
    LO k = v / nvxy;
    LO i = ij % nvx;
    LO j = ij / nvx;
    coords[v * 3 + 0] = i * dx;
    coords[v * 3 + 1] = j * dy;
    coords[v * 3 + 2] = k * dz;
  };
  parallel_for(nv, fill_coords);
  Write<LO> hv2v(nh * 8);
  auto fill_conn = LAMBDA(LO h) {
    LO ij = h % nxy;
    LO k = h / nxy;
    LO i = ij % nx;
    LO j = ij / nx;
    hv2v[h * 8 + 0] = (k + 0) * nvxy + (j + 0) * nvx + (i + 0);
    hv2v[h * 8 + 1] = (k + 0) * nvxy + (j + 0) * nvx + (i + 1);
    hv2v[h * 8 + 2] = (k + 0) * nvxy + (j + 1) * nvx + (i + 1);
    hv2v[h * 8 + 3] = (k + 0) * nvxy + (j + 1) * nvx + (i + 0);
    hv2v[h * 8 + 4] = (k + 1) * nvxy + (j + 0) * nvx + (i + 0);
    hv2v[h * 8 + 5] = (k + 1) * nvxy + (j + 0) * nvx + (i + 1);
    hv2v[h * 8 + 6] = (k + 1) * nvxy + (j + 1) * nvx + (i + 1);
    hv2v[h * 8 + 7] = (k + 1) * nvxy + (j + 1) * nvx + (i + 0);
  };
  parallel_for(nh, fill_conn);
  *hv2v_out = hv2v;
  *coords_out = coords;
}

template <Int dim>
static Read<I32> box_centroids_class_ids(Reals centroids,
    Few<LO, 3> nel, Vector<3> l) {
  CHECK(centroids.size() % dim == 0);
  auto npts = centroids.size() / dim;
  Vector<dim> dists;
  for (Int i = 0; i < dim; ++i) dists[i] = l[i] / (nel[i] * 8);
  auto class_ids = Write<I32>(npts);
  auto f = LAMBDA(Int i) {
    auto x = get_vector<dim>(centroids, i);
    Int id = 0;
    for (Int j = dim - 1; j >= 0; --j) {
      id *= 3;
      if (x[j] > (l[j] - dists[j])) id += 2;
      else if (x[j] > dists[j]) id += 1;
    }
    class_ids[i] = id;
  };
  parallel_for(npts, f);
  return class_ids;
}

void set_box_class_ids(Mesh* mesh,
    Real x, Real y, Real z,
    LO nx, LO ny, LO nz) {
  Few<LO, 3> nel({nx, ny, nz});
  Vector<3> l({x, y, z});
  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    Reals centroids;
    if (ent_dim) {
      centroids = average_field(mesh, ent_dim,
          LOs(mesh->nents(ent_dim), 0, 1), mesh->dim(), mesh->coords());
    } else {
      centroids = mesh->coords();
    }
    Read<I32> class_ids;
    if (mesh->dim() == 3)
      class_ids = box_centroids_class_ids<3>(centroids, nel, l);
    if (mesh->dim() == 2)
      class_ids = box_centroids_class_ids<2>(centroids, nel, l);
    mesh->add_tag<I32>(ent_dim, "class_id", 1, OMEGA_H_INHERIT,
        OMEGA_H_DO_OUTPUT, class_ids);
  }
}

}  // end namespace Omega_h
