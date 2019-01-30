#ifndef OMEGA_H_TRANSFER_FACE_HPP
#define OMEGA_H_TRANSFER_FACE_HPP

#include <Omega_h_adapt.hpp>
#include <Omega_h_adj.hpp>
#include <Omega_h_tag.hpp>

namespace Omega_h {
  template <typename T>
  void transer_div_free_face_flux(Mesh *src_mesh,
				  Mesh *target_mesh,
				  Read<T> src_data,
				  Read<T> target_data);

#define INST_DECL(T)                                                    \
  extern template void transer_div_free_face_flux(Mesh *src_mesh,       \
						  Mesh *target_mesh,	\
						  Read<T> src_data,	\
						  Read<T> targeta_data);	\
INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL

}  // end namespace Omega_h

#endif
