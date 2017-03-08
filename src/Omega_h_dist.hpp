#ifndef OMEGA_H_DIST_HPP
#define OMEGA_H_DIST_HPP

#include <Omega_h_comm.hpp>
#include <Omega_h_remotes.hpp>

namespace Omega_h {

class Dist {
  CommPtr parent_comm_;
  LOs roots2items_[2];
  LOs items2content_[2];
  LOs msgs2content_[2];
  CommPtr comm_[2];

 public:
  Dist();
  Dist(Dist const& other);
  Dist& operator=(Dist const& other);
  Dist(CommPtr comm, Remotes fitems2rroots, LO nrroots);
  void set_parent_comm(CommPtr parent_comm);
  void set_dest_ranks(Read<I32> items2ranks);
  void set_dest_idxs(LOs fitems2rroots, LO nrroots);
  void set_roots2items(LOs froots2fitems);
  Dist invert() const;
  template <typename T>
  Read<T> exch(Read<T> data, Int width) const;
  template <typename T>
  Read<T> exch_reduce(Read<T> data, Int width, Omega_h_Op op) const;
  CommPtr parent_comm() const;
  CommPtr comm() const;
  LOs msgs2content() const;
  LOs content2msgs() const;
  LOs items2msgs() const;
  LOs roots2items() const;
  Read<I32> msgs2ranks() const;
  Read<I32> items2ranks() const;
  LOs items2dest_idxs() const;
  Remotes items2dests() const;
  LO nitems() const;
  LO nroots() const;
  LO ndests() const;
  LO nsrcs() const;
  void change_comm(CommPtr new_comm);
  Remotes exch(Remotes data, Int width) const;

 private:
  void copy(Dist const& other);
  enum { F, R };
};

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template Read<T> Dist::exch(Read<T> data, Int width) const;           \
  extern template Read<T> Dist::exch_reduce<T>(                                \
      Read<T> data, Int width, Omega_h_Op op) const;
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

}

#endif
