#ifndef OMEGA_H_HPP
#define OMEGA_H_HPP

#include <iosfwd>
#include <new>
#include <vector>

#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_adapt.hpp>

namespace Omega_h {

void build_from_elems2verts(
    Mesh* mesh, CommPtr comm, Int edim, LOs ev2v, Read<GO> vert_globals);
void build_from_elems2verts(Mesh* mesh, Int edim, LOs ev2v, LO nverts);
void build_from_elems_and_coords(Mesh* mesh, Int edim, LOs ev2v, Reals coords);
void build_box(Mesh* mesh, Real x, Real y, Real z, LO nx, LO ny, LO nz);

void classify_by_angles(Mesh* mesh, Real sharp_angle);

Real repro_sum(Reals a);
Real repro_sum(CommPtr comm, Reals a);
void repro_sum(CommPtr comm, Reals a, Int ncomps, Real result[]);
Real repro_sum_owned(Mesh* mesh, Int dim, Reals a);

OMEGA_H_INLINE bool code_is_flipped(I8 code) { return code & 1; }

OMEGA_H_INLINE Int code_rotation(I8 code) { return (code >> 1) & 3; }

OMEGA_H_INLINE Int code_which_down(I8 code) { return (code >> 3); }

Read<I8> mark_class_closure(
    Mesh* mesh, Int ent_dim, Int class_dim, I32 class_id);
Read<I8> mark_class_closures(Mesh* mesh, Int ent_dim,
    std::vector<Int> class_dims, std::vector<I32> class_ids);
void fix_momentum_velocity_verts(
    Mesh* mesh, Int class_dim, I32 class_id, Int comp);

bool warp_to_limit(Mesh* mesh, AdaptOpts const& opts);
bool approach_size_field(Mesh* mesh, AdaptOpts const& opts);

Reals find_implied_size(Mesh* mesh);
Reals find_implied_metric(Mesh* mesh);
void axes_from_metric_field(
    Mesh* mesh, std::string const& metric_name, std::string const& axis_prefix);
Reals limit_size_field_gradation(
    Mesh* mesh, Reals values, Real max_rate, Real tol = 1e-3);
Reals expected_elems_per_elem_iso(Mesh* mesh, Reals v2h);
Reals expected_elems_per_elem_metric(Mesh* mesh, Reals v2m);
Real size_scalar_for_nelems(Mesh* mesh, Reals v2h, Real target_nelems);
Real metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems);
Reals smooth_metric_once(Mesh* mesh, Reals v2m);
Reals smooth_isos_once(Mesh* mesh, Reals v2h);
Reals get_curvature_isos(Mesh* mesh, Real segment_angle, Real max_size);
Reals get_gradient_isos(
    Mesh* mesh, Real error_bound, Real max_size, Reals scalar_field);
Reals clamp_deforming_isos(Mesh* mesh, Reals isos, Real min_size,
    Real max_interior_size, Real max_boundary_size);

Reals recover_hessians(Mesh* mesh, Reals vert_values);
Reals metric_from_hessians(Int dim, Reals hessians, Real eps, Real hmax);
Reals metric_for_nelems_from_hessians(
    Mesh* mesh, Real target_nelems, Real tolerance, Reals hessians, Real hmax);

template <typename T, Int n>
class Few {
  using UninitT = typename std::aligned_storage<sizeof(T), alignof(T)>::type;
  UninitT array_[n];

 public:
  enum { size = n };
  OMEGA_H_INLINE T* data() { return reinterpret_cast<T*>(array_); }
  OMEGA_H_INLINE T const* data() const {
    return reinterpret_cast<T const*>(array_);
  }
  OMEGA_H_INLINE T volatile* data() volatile {
    return reinterpret_cast<T volatile*>(array_);
  }
  OMEGA_H_INLINE T const volatile* data() const volatile {
    return reinterpret_cast<T const volatile*>(array_);
  }
  OMEGA_H_INLINE T& operator[](Int i) { return data()[i]; }
  OMEGA_H_INLINE T const& operator[](Int i) const { return data()[i]; }
  OMEGA_H_INLINE T volatile& operator[](Int i) volatile { return data()[i]; }
  OMEGA_H_INLINE T const volatile& operator[](Int i) const volatile {
    return data()[i];
  }
  Few(std::initializer_list<T> l) {
    Int i = 0;
    for (auto it = l.begin(); it != l.end(); ++it) {
      new (data() + (i++)) T(*it);
    }
  }
  OMEGA_H_INLINE Few() {
    for (Int i = 0; i < n; ++i) new (data() + i) T();
  }
  OMEGA_H_INLINE ~Few() {
    for (Int i = 0; i < n; ++i) (data()[i]).~T();
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const& rhs) volatile {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE Few(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
  OMEGA_H_INLINE Few(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
};

template <typename T>
OMEGA_H_INLINE T max2(T a, T b) {
  return (b > a) ? (b) : (a);
}

template <typename T>
OMEGA_H_INLINE T min2(T a, T b) {
  return (b < a) ? (b) : (a);
}

template <typename T>
OMEGA_H_INLINE void swap2(T& a, T& b) {
  T c = a;
  a = b;
  b = c;
}

bool ends_with(std::string const& s, std::string const& suffix);

}  // end namespace Omega_h

#endif
