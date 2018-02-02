#ifndef OMEGA_H_OVERLAY_HPP
#define OMEGA_H_OVERLAY_HPP

#include <Omega_h_mesh.hpp>

namespace Omega_h {

class Overlay {
 public:
  Overlay(int* argc, char*** argv, const Omega_h::Vector<3>& overlayGridCenter,
      double cellSize, size_t Nx, size_t Ny, size_t Nz
#ifdef OMEGA_H_USE_MPI
      ,
      MPI_Comm comm
#endif
  );

  std::array<size_t, 2> get_face_cells(size_t face) const;

  const std::vector<size_t>& get_cells() const;

  const std::vector<size_t>& get_edges() const;

  const std::vector<size_t>& get_faces() const;

  const std::vector<size_t>& get_nodes() const;

  std::array<size_t, 8> get_cell_nodes(size_t cell) const;

  std::array<size_t, 6> get_cell_faces(size_t cell) const;

  std::array<size_t, 4> get_face_nodes(size_t face) const;

  std::array<size_t, 2> get_edge_nodes(size_t edge) const;

  Omega_h::Vector<3> get_node_coordinates(size_t node) const;

  Omega_h::Vector<3> get_cell_center_location(size_t cell) const;

  size_t get_invalid_cell_handle() const;

  double get_cell_size() const;

 private:
  double cell_size;
  Omega_h::Library library;
  Omega_h::Mesh mesh;
  std::vector<size_t> silly_cells;
  std::vector<size_t> silly_edges;
  std::vector<size_t> silly_faces;
  std::vector<size_t> silly_nodes;
};

}  // namespace Omega_h

#endif
