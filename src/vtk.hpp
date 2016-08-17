#ifndef VTK_HPP
#define VTK_HPP

#include "internal.hpp"

namespace Omega_h {

namespace vtk {

std::string get_pvtu_path(std::string const& step_path);

std::string get_pvd_path(std::string const& root_path);

void write_pvtu(std::ostream& stream, Mesh* mesh, Int cell_dim,
    std::string const& piecepath);

void write_pvtu(std::string const& filename, Mesh* mesh, Int cell_dim,
    std::string const& piecepath);

std::streampos write_initial_pvd(std::string const& root_path);
void update_pvd(std::string const& root_path, std::streampos* pos_inout,
    Int step, Real time);

void read_vtu(std::istream& stream, CommPtr comm, Mesh* mesh);

void read_pvtu(std::istream& stream, CommPtr comm, I32* npieces_out,
    std::string* vtupath_out);

void read_pvtu(std::string const& pvtupath, CommPtr comm, I32* npieces_out,
    std::string* vtupath_out);

void read_parallel(std::string const& pvtupath, CommPtr comm, Mesh* mesh);

void read_pvd(std::istream& stream, std::vector<Real>* times_out,
    std::vector<std::string>* pvtupaths_out);

void read_pvd(std::string const& pvdpath, std::vector<Real>* times_out,
    std::vector<std::string>* pvtupaths_out);
}

}  // end namespace Omega_h

#endif
