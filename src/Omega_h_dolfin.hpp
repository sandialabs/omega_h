#ifndef OMEGA_H_DOLFIN_HPP
#define OMEGA_H_DOLFIN_HPP

#include <Omega_h_config.h>

#ifdef OMEGA_H_USE_MPI
#define HAS_MPI // omg DOLFIN this is so bad
#endif
#include <dolfin.h>
#ifdef OMEGA_H_USE_MPI
#undef HAS_MPI // omg DOLFIN this is so bad
#endif

namespace Omega_h {

class Mesh;

void to_dolfin(dolfin::Mesh& mesh_dolfin, Mesh* mesh_osh);
void from_dolfin(Mesh* mesh_osh, dolfin::Mesh const& mesh_dolfin);

}

#endif // OMEGA_H_DOLFIN_HPP
