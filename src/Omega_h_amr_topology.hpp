#ifndef OMEGA_H_AMR_TOPOLOGY_HPP
#define OMEGA_H_AMR_TOPOLOGY_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_few.hpp>
#include <Omega_h_template_up.hpp>

namespace Omega_h {

class Mesh;

void mark_amr(Mesh* mesh, Read<Byte> elem_mark);

Few<LO, 4> count_amr(Mesh* mesh);

}

#endif
