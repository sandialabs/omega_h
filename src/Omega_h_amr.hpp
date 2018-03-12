#ifndef OMEGA_H_AMR_HPP
#define OMEGA_H_AMR_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

void amr_refine(Mesh* mesh, Bytes elems_are_marked);

}

#endif
