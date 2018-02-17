#ifndef OMEGA_H_DOLFIN_HPP
#define OMEGA_H_DOLFIN_HPP

#include <Omega_h_config.h>
#include "Omega_h_mesh.hpp"

#if defined(OMEGA_H_USE_MPI) && (!defined(HAS_MPI))
#define HAS_MPI
#define OMEGA_H_SET_HAS_MPI
#endif

#include <dolfin.h>

#ifdef OMEGA_H_SET_HAS_MPI
#undef HAS_MPI
#endif

namespace Omega_h {

class Mesh;

void to_dolfin(dolfin::Mesh& mesh_dolfin, Mesh* mesh_osh);
void from_dolfin(Mesh* mesh_osh, dolfin::Mesh const& mesh_dolfin);
void from_dolfin(
    Mesh* mesh_osh, dolfin::Function const& function, std::string const& name);

}  // namespace Omega_h

#endif  // OMEGA_H_DOLFIN_HPP
