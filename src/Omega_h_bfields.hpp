#ifndef OMEGA_H_BFIELDS_HPP
#define OMEGA_H_BFIELDS_HPP

#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <Omega_h_adj.hpp>

namespace Omega_h {

class Mesh {
 public:

  Adj ask_revClass (Int edim);
  Adj ask_revClass (Int edim, LOs g_ids);
  bool has_revClass (Int edim) const;
  Adj get_revClass (Int edim) const;
  Adj ask_revClass_downAdj (Int from, Int to);
  Adj derive_revClass (Int edim);

  template <typename T>
  void add_boundaryField(Int ent_dim, std::string const& name, Int ncomps);
  template <typename T>
  void add_boundaryField(
    Int ent_dim, std::string const& name, Int ncomps, Read<T> array,
    bool internal = false);
  template <typename T>
  Read<T> get_boundaryField_array(Int ent_dim, std::string const& name) const;
  template <typename T>
  void set_boundaryField_array(
    Int ent_dim, std::string const& name, Read<T> array, bool internal = false);
  void reduce_boundaryField(Int ent_dim, std::string const& name,
     Omega_h_Op op);
  void sync_boundaryField(Int ent_dim, std::string const& name);
  bool has_boundaryField(Int ent_dim, std::string const& name) const;
  void remove_boundaryField(Int ent_dim, std::string const& name);

  template <typename T>
  void change_tagToBoundary(Int ent_dim, Int ncomps, std::string const& name);
  template <typename T>
  void change_tagToMesh(Int ent_dim, Int ncomps, std::string const& name);

  void change_all_bFieldsToMesh();
  void change_all_bFieldsToBoundary();
  bool has_anyBoundaryField();
  bool has_allMeshTags();

 private:

  AdjPtr revClass_[DIMS];

};

#ifdef OMEGA_H_USE_CUDA
__host__
#endif
    void
    assign(Mesh& a, Mesh const& b);

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template void Mesh::change_tagToBoundary<T>(                          \
      Int ent_dim, Int ncomps, std::string const& name);                       \
  extern template void Mesh::change_tagToMesh<T>(                              \
      Int ent_dim, Int ncomps, std::string const& name);                       \
  extern template Read<T> Mesh::get_boundaryField_array<T>(                    \
      Int dim, std::string const& name) const;                                 \
  extern template void Mesh::add_boundaryField<T>(                             \
      Int dim, std::string const& name, Int ncomps);                           \
  extern template void Mesh::add_boundaryField<T>(                             \
      Int dim, std::string const& name, Int ncomps, Read<T> array,             \
      bool internal);                                                          \
  extern template void Mesh::set_boundaryField_array(                          \
      Int dim, std::string const& name, Read<T> array, bool internal);         
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

}  // namespace Omega_h

#endif
