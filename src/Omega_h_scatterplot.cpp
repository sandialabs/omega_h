#include <Omega_h_scatterplot.hpp>

#include <Omega_h_for.hpp>
#include <Omega_h_linpart.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_owners.hpp>

#include <fstream>
#include <iomanip>

namespace Omega_h {

template <Int dim>
Reals get_linear_scatter_coords(
    Reals coords, Vector<dim> direction, Vector<dim> origin) {
  auto offset = origin * direction;
  auto n = divide_no_remainder(coords.size(), dim);
  auto coords_1d_w = Write<Real>(n);
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto x = get_vector<dim>(coords, i);
    auto d = x * direction - offset;
    coords_1d_w[i] = d;
  };
  parallel_for(n, f);
  return coords_1d_w;
}

template <Int dim>
Reals get_radial_scatter_coords(Reals coords, Vector<dim> center) {
  auto n = divide_no_remainder(coords.size(), dim);
  auto coords_1d_w = Write<Real>(n);
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto x = get_vector<dim>(coords, i);
    auto d = norm(x - center);
    coords_1d_w[i] = d;
  };
  parallel_for(n, f);
  return coords_1d_w;
}

void write_scatterplot(std::string const& path, CommPtr comm, Reals coords_1d,
    Reals data, std::string const& separator) {
  HostRead<Real> coords_1d_h(coords_1d);
  HostRead<Real> data_h(data);
  auto ncomps = divide_no_remainder(data.size(), coords_1d.size());
  for (I32 rank = 0; rank < comm->size(); ++rank) {
    std::ofstream os;
    if (rank)
      os.open(path.c_str(), std::ios_base::ate);
    else
      os.open(path.c_str(), std::ios_base::trunc);
    os << std::scientific << std::setprecision(17);
    if (comm->rank() == rank) {
      for (LO i = 0; i < coords_1d_h.size(); ++i) {
        os << coords_1d_h[i];
        for (Int j = 0; j < ncomps; ++j) {
          os << separator << data_h[i * ncomps + j];
        }
        os << '\n';
      }
    }
    os.close();
    comm->barrier();
  }
}

void write_scatterplot(std::string const& path, Mesh* mesh, Int ent_dim,
    Reals coords_1d, Reals data, std::string const& separator) {
  auto comm = mesh->comm();
  auto ncomps = divide_no_remainder(data.size(), coords_1d.size());
  if (comm->size() > 1) {
    auto copies2lins = copies_to_linear_owners(comm, mesh->globals(ent_dim));
    coords_1d = reduce_data_to_owners(coords_1d, copies2lins, 1);
    data = reduce_data_to_owners(data, copies2lins, ncomps);
  }
  write_scatterplot(path, comm, coords_1d, data, separator);
}

template <Int dim>
void write_linear_scatterplot(std::string const& path, Mesh* mesh, Int ent_dim,
    Reals data, Vector<dim> direction, Vector<dim> origin,
    std::string const& separator) {
  auto coords =
      average_field(mesh, ent_dim, dim, mesh->coords());  // get centroids
  auto coords_1d = get_linear_scatter_coords(coords, direction, origin);
  write_scatterplot(path, mesh, ent_dim, coords_1d, data, separator);
}

template <Int dim>
void write_radial_scatterplot(std::string const& path, Mesh* mesh, Int ent_dim,
    Reals data, Vector<dim> center, std::string const& separator) {
  auto coords =
      average_field(mesh, ent_dim, dim, mesh->coords());  // get centroids
  auto coords_1d = get_radial_scatter_coords(coords, center);
  write_scatterplot(path, mesh, ent_dim, coords_1d, data, separator);
}

#define OMEGA_H_EXPL_INST(dim)                                                 \
  template Reals get_linear_scatter_coords(                                    \
      Reals coords, Vector<dim> direction, Vector<dim> origin);                \
  template Reals get_radial_scatter_coords(Reals coords, Vector<dim> center);  \
  template void write_linear_scatterplot(std::string const& path, Mesh* mesh,  \
      Int ent_dim, Reals data, Vector<dim> direction, Vector<dim> origin,      \
      std::string const& separator);                                           \
  template void write_radial_scatterplot(std::string const& path, Mesh* mesh,  \
      Int ent_dim, Reals data, Vector<dim> center,                             \
      std::string const& separator);
OMEGA_H_EXPL_INST(1)
OMEGA_H_EXPL_INST(2)
OMEGA_H_EXPL_INST(3)
#undef OMEGA_H_EXPL_INST

}  // namespace Omega_h
