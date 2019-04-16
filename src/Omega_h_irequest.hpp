#ifndef OMEGA_H_  ASYNC_REQUEST_HPP
#define OMEGA_H_  ASYNC_REQUEST_HPP

#include "mpi.h"

namespace Omega_h {

template <typename T>
class IRequest {
 public:
  IRequest(MPI_request* request, const std::function<Read<T>()>& future)
 : request_(request), future(future) {}

 bool completed();
 Read<T> get();

 private:
  MPI_Request* request;
  const std::function<T()>& future;
};

#define OMEGA_H_INST(T)                                                        \
  template Tag<T> const* Mesh::get_tag<T>(Int dim, std::string const& name)    \
      const;                                                                   \
  template Read<T> Mesh::get_array<T>(Int dim, std::string const& name) const; \
  template void Mesh::add_tag<T>(                                              \
      Int dim, std::string const& name, Int ncomps);                           \
  template void Mesh::add_tag<T>(Int dim, std::string const& name, Int ncomps, \
      Read<T> array, bool internal);                                           \
  template void Mesh::set_tag(                                                 \
      Int dim, std::string const& name, Read<T> array, bool internal);         \
  template Read<T> Mesh::sync_array(Int ent_dim, Read<T> a, Int width);        \
  template Read<T> Mesh::owned_array(Int ent_dim, Read<T> a, Int width);       \
  template Read<T> Mesh::sync_subset_array(                                    \
      Int ent_dim, Read<T> a_data, LOs a2e, T default_val, Int width);         \
  template Read<T> Mesh::reduce_array(                                         \
      Int ent_dim, Read<T> a, Int width, Omega_h_Op op);
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST


} // namespace Omega_h

#endif //OMEGA_H_  ASYNC_REQUEST_HPP
