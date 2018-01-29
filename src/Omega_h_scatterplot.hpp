#ifndef OMEGA_H_SCATTERPLOT_HPP
#define OMEGA_H_SCATTERPLOT_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_vector.hpp>

#include <string>

namespace Omega_h {

class Mesh;

template <Int dim>
Reals get_linear_scatter_coords(
    Reals coords, Vector<dim> direction, Vector<dim> origin);
template <Int dim>
Reals get_radial_scatter_coords(Reals coords, Vector<dim> center);
void write_scatterplot(std::string const& path, CommPtr comm, Reals coords_1d,
    Reals data, std::string const& separator);
void write_scatterplot(std::string const& path, Mesh* mesh, Int ent_dim,
    Reals coords_1d, Reals data, std::string const& separator);
template <Int dim>
void write_linear_scatterplot(std::string const& path, Mesh* mesh, Int ent_dim,
    Reals data, Vector<dim> direction, Vector<dim> origin,
    std::string const& separator);
template <Int dim>
void write_radial_scatterplot(std::string const& path, Mesh* mesh, Int ent_dim,
    Reals data, Vector<dim> center, std::string const& separator);

#define OMEGA_H_EXPL_INST_DECL(dim)                                            \
  extern template Reals get_linear_scatter_coords(                             \
      Reals coords, Vector<dim> direction, Vector<dim> origin);                \
  extern template Reals get_radial_scatter_coords(                             \
      Reals coords, Vector<dim> center);                                       \
  extern template void write_linear_scatterplot(std::string const& path,       \
      Mesh* mesh, Int ent_dim, Reals data, Vector<dim> direction,              \
      Vector<dim> origin, std::string const& separator);                       \
  extern template void write_radial_scatterplot(std::string const& path,       \
      Mesh* mesh, Int ent_dim, Reals data, Vector<dim> center,                 \
      std::string const& separator);
OMEGA_H_EXPL_INST_DECL(1)
OMEGA_H_EXPL_INST_DECL(2)
OMEGA_H_EXPL_INST_DECL(3)
#undef OMEGA_H_EXPL_INST_DECL

}  // end namespace Omega_h

#endif
