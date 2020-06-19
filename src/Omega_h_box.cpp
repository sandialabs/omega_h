#include "Omega_h_box.hpp"

#include "Omega_h_few.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_vector.hpp"

namespace Omega_h {

void make_1d_box(Real x, LO nx, LOs* ev2v_out, Reals* coords_out) {
  LO ne = nx;
  LO nv = nx + 1;
  Real dx = x / nx;
  Write<Real> coords(nv);
  auto fill_coords = OMEGA_H_LAMBDA(LO v) {
    LO i = v % nv;
    coords[v] = i * dx;
  };
  parallel_for(nv, fill_coords, "make_1d_box(coords)");
  Write<LO> ev2v(ne * 2);
  auto fill_conn = OMEGA_H_LAMBDA(LO q) {
    LO i = q % ne;
    ev2v[q * 2 + 0] = i + 0;
    ev2v[q * 2 + 1] = i + 1;
  };
  parallel_for(ne, fill_conn, "make_1d_box(conn)");
  *ev2v_out = ev2v;
  *coords_out = coords;
}

void make_2d_box(
    Real x, Real y, LO nx, LO ny, LOs* qv2v_out, Reals* coords_out) {
  LO nq = nx * ny;
  LO nvx = nx + 1;
  LO nvy = ny + 1;
  LO nv = nvx * nvy;
  Real dx = x / nx;
  Real dy = y / ny;
  Write<Real> coords(nv * 2);
  auto fill_coords = OMEGA_H_LAMBDA(LO v) {
    LO i = v % nvx;
    LO j = v / nvx;
    coords[v * 2 + 0] = i * dx;
    coords[v * 2 + 1] = j * dy;
  };
  parallel_for(nv, fill_coords, "make_2d_box(coords)");
  Write<LO> qv2v(nq * 4);
  auto fill_conn = OMEGA_H_LAMBDA(LO q) {
    LO i = q % nx;
    LO j = q / nx;
    qv2v[q * 4 + 0] = (j + 0) * nvx + (i + 0);
    qv2v[q * 4 + 1] = (j + 0) * nvx + (i + 1);
    qv2v[q * 4 + 2] = (j + 1) * nvx + (i + 1);
    qv2v[q * 4 + 3] = (j + 1) * nvx + (i + 0);
  };
  parallel_for(nq, fill_conn, "make_2d_box(conn)");
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
  auto fill_coords = OMEGA_H_LAMBDA(LO v) {
    LO ij = v % nvxy;
    LO k = v / nvxy;
    LO i = ij % nvx;
    LO j = ij / nvx;
    coords[v * 3 + 0] = i * dx;
    coords[v * 3 + 1] = j * dy;
    coords[v * 3 + 2] = k * dz;
  };
  parallel_for(nv, fill_coords, "make_3d_box(coords)");
  Write<LO> hv2v(nh * 8);
  auto fill_conn = OMEGA_H_LAMBDA(LO h) {
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
  parallel_for(nh, fill_conn, "make_3d_box(conn)");
  *hv2v_out = hv2v;
  *coords_out = coords;
}

template <Int dim>
void classify_box_dim(Mesh* mesh, Int ent_dim, Reals centroids,
    Few<LO, 3> nel, Vector<3> const l) {
  OMEGA_H_CHECK(centroids.size() % dim == 0);
  auto npts = centroids.size() / dim;
  Vector<dim> dists;
  /* we assume that if an entity should not be classified on
     the boundary surface, its centroid is more than an (1/32)
     of a cell width away from said boundary */
  for (Int i = 0; i < dim; ++i) dists[i] = l[i] / (nel[i] * 32);
  auto class_ids = Write<ClassId>(npts);
  auto class_dims = Write<Byte>(npts);
  auto f = OMEGA_H_LAMBDA(Int i) {
    auto x = get_vector<dim>(centroids, i);
    Int id = 0;
    Int class_dim = 0;
    for (Int j = dim - 1; j >= 0; --j) {
      id *= 3;
      if (x[j] > (l[j] - dists[j])) {
        /* case 1: point lies on the upper boundary */
        id += 2;
      } else if (x[j] > dists[j]) {
        /* case 2: point lies on the interior */
        id += 1;
        ++class_dim;
      }
    }
    class_ids[i] = id;
    class_dims[i] = Byte(class_dim);
  };
  parallel_for(npts, f, "set_box_class_ids");
  mesh->add_tag<ClassId>(ent_dim, "class_id", 1, class_ids);
  mesh->add_tag<Byte>(ent_dim, "class_dim", 1, class_dims);
}

void classify_box(Mesh* mesh, Real x, Real y, Real z, LO nx, LO ny, LO nz) {
  Few<LO, 3> nel({nx, ny, nz});
  Vector<3> l({x, y, z});
  for (Int ent_dim = 0; ent_dim <= mesh->dim(); ++ent_dim) {
    Reals centroids;
    if (ent_dim) {
      centroids = average_field(mesh, ent_dim, LOs(mesh->nents(ent_dim), 0, 1),
          mesh->dim(), mesh->coords());
    } else {
      centroids = mesh->coords();
    }
    Read<LO> class_ids;
    if (mesh->dim() == 3)
      classify_box_dim<3>(mesh, ent_dim, centroids, nel, l);
    else if (mesh->dim() == 2)
      classify_box_dim<2>(mesh, ent_dim, centroids, nel, l);
    else if (mesh->dim() == 1)
      classify_box_dim<1>(mesh, ent_dim, centroids, nel, l);
    else
      Omega_h_fail("classify_box: dimension isn't 1, 2, or 3!");
  }
}

ClassSets get_box_class_sets(Int dim) {
  ClassSets sets;
  if (dim == 1) {
    sets["x-"] = {{0, 0}};
    sets["body"] = {{1, 1}};
    sets["x+"] = {{0, 2}};
  } else if (dim == 2) {
    sets["y-"] = {{1, 1}};
    sets["x-"] = {{1, 3}};
    sets["body"] = {{2, 4}};
    sets["x+"] = {{1, 5}};
    sets["y+"] = {{1, 7}};
  } else if (dim == 3) {
    sets["z-"] = {{2, 4}};
    sets["y-"] = {{2, 10}};
    sets["x-"] = {{2, 12}};
    sets["body"] = {{3, 13}};
    sets["x+"] = {{2, 14}};
    sets["y+"] = {{2, 16}};
    sets["z+"] = {{2, 22}};
  }
  return sets;
}

}  // end namespace Omega_h
