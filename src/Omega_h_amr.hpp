#ifndef OMEGA_H_AMR_HPP
#define OMEGA_H_AMR_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_adapt.hpp>

namespace Omega_h {

class Mesh;

void amr_refine(Mesh* mesh, Bytes elems_are_marked, TransferOpts xfer_opts);

}

#endif
