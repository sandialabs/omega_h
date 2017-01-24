#ifndef SIMPLICES_HPP
#define SIMPLICES_HPP

#include "internal.hpp"

namespace Omega_h {

INLINE Int down_template(Int elem_dim, Int bdry_dim, Int which_bdry, Int which_vert) {
  switch (elem_dim) {
    case 2:
      switch (bdry_dim) {
        case 0: return which_bdry;
        case 1:
          switch (which_bdry) {
            case 0:
              switch (which_vert) {
                case 0: return 0;
                case 1: return 1;
              }
            case 1:
              switch (which_vert) {
                case 0: return 1;
                case 1: return 2;
              }
            case 2:
              switch (which_vert) {
                case 0: return 2;
                case 1: return 0;
              }
          }
      }
    case 3:
      switch (bdry_dim) {
        case 0: return which_bdry;
        case 1:
          switch (which_bdry) {
            case 0:
              switch (which_vert) {
                case 0: return 0;
                case 1: return 1;
              }
            case 1:
              switch (which_vert) {
                case 0: return 1;
                case 1: return 2;
              }
            case 2:
              switch (which_vert) {
                case 0: return 2;
                case 1: return 0;
              }
            case 3:
              switch (which_vert) {
                case 0: return 0;
                case 1: return 3;
              }
            case 4:
              switch (which_vert) {
                case 0: return 1;
                case 1: return 3;
              }
            case 5:
              switch (which_vert) {
                case 0: return 2;
                case 1: return 3;
              }
          }
        case 2:
          switch (which_bdry) {
            case 0:
              switch (which_vert) {
                case 0: return 0;
                case 1: return 2;
                case 2: return 1;
              }
            case 1:
              switch (which_vert) {
                case 0: return 0;
                case 1: return 1;
                case 2: return 3;
              }
            case 2:
              switch (which_vert) {
                case 0: return 1;
                case 1: return 2;
                case 2: return 3;
              }
            case 3:
              switch (which_vert) {
                case 0: return 2;
                case 1: return 0;
                case 2: return 3;
              }
          }
      }
  }
  return -1;
}

struct TemplateUp {
  Int up;
  Int which_down;
  bool is_flipped;
};
CONSTANT static TemplateUp const fve0[] = {{0, 0, 0}, {2, 1, 0}};
CONSTANT static TemplateUp const fve1[] = {{1, 0, 0}, {0, 1, 0}};
CONSTANT static TemplateUp const fve2[] = {{2, 0, 0}, {1, 1, 0}};
CONSTANT static TemplateUp const* const fve_[] = {fve0, fve1, fve2};
CONSTANT static TemplateUp const* const* const f_u_[] = {fve_};
CONSTANT static TemplateUp const rve0[] = {{0, 0, 0}, {2, 1, 0}, {3, 0, 0}};
CONSTANT static TemplateUp const rve1[] = {{1, 0, 0}, {0, 1, 0}, {4, 0, 0}};
CONSTANT static TemplateUp const rve2[] = {{2, 0, 0}, {1, 1, 0}, {5, 0, 0}};
CONSTANT static TemplateUp const rve3[] = {{5, 1, 0}, {4, 1, 0}, {3, 1, 0}};
CONSTANT static TemplateUp const* const rve_[] = {rve0, rve1, rve2, rve3};
CONSTANT static TemplateUp const ref0[] = {{0, 2, 1}, {1, 0, 0}};
CONSTANT static TemplateUp const ref1[] = {{0, 1, 1}, {2, 0, 0}};
CONSTANT static TemplateUp const ref2[] = {{0, 0, 1}, {3, 0, 0}};
CONSTANT static TemplateUp const ref3[] = {{1, 2, 1}, {3, 1, 0}};
CONSTANT static TemplateUp const ref4[] = {{2, 2, 1}, {1, 1, 0}};
CONSTANT static TemplateUp const ref5[] = {{3, 2, 1}, {2, 1, 0}};
CONSTANT static TemplateUp const* const ref_[] = {
    ref0, ref1, ref2, ref3, ref4, ref5};
CONSTANT static TemplateUp const* const* const r_u_[] = {rve_, ref_};
CONSTANT static TemplateUp const* const* const* const up_templates[] = {
    0, 0, f_u_, r_u_};

/* workaround a compiler bug in CUDA, see DownTemplate<> */
template <Int hdim, Int ldim>
struct UpTemplate;
template <>
struct UpTemplate<3, 1> {
  DEVICE static TemplateUp get(Int a, Int b) { return ref_[a][b]; }
};
template <>
struct UpTemplate<3, 0> {
  DEVICE static TemplateUp get(Int a, Int b) { return rve_[a][b]; }
};
template <>
struct UpTemplate<2, 0> {
  DEVICE static TemplateUp get(Int a, Int b) { return fve_[a][b]; }
};

CONSTANT static Int const feov[] = {1, 2, 0};
CONSTANT static Int const fvoe[] = {2, 0, 1};
CONSTANT static Int const* const fo_[] = {feov, fvoe};
CONSTANT static Int const rfov[] = {2, 3, 1, 0};
CONSTANT static Int const rvof[] = {3, 2, 0, 1};
CONSTANT static Int const reoe[] = {5, 3, 4, 1, 2, 0};
CONSTANT static Int const* const ro_[] = {rfov, reoe, rvof};
CONSTANT static Int const* const* const opposite_templates[] = {0, 0, fo_, ro_};

/* workaround a compiler bug in CUDA, see DownTemplate<> */
template <Int hdim, Int ldim>
struct OppositeTemplate;
template <>
struct OppositeTemplate<3, 2> {
  DEVICE static Int get(Int a) { return rvof[a]; }
};
template <>
struct OppositeTemplate<3, 1> {
  DEVICE static Int get(Int a) { return reoe[a]; }
};
template <>
struct OppositeTemplate<3, 0> {
  DEVICE static Int get(Int a) { return rfov[a]; }
};
template <>
struct OppositeTemplate<2, 1> {
  DEVICE static Int get(Int a) { return fvoe[a]; }
};
template <>
struct OppositeTemplate<2, 0> {
  DEVICE static Int get(Int a) { return feov[a]; }
};

extern Int const simplex_degrees[DIMS][DIMS];
extern char const* const singular_names[DIMS];
extern char const* const plural_names[DIMS];

template <Int dim, Int low, Int high>
struct AvgDegree;

template <>
struct AvgDegree<2, 0, 1> {
  static constexpr Int value = 6;
};

template <>
struct AvgDegree<2, 0, 2> {
  static constexpr Int value = 6;
};

template <>
struct AvgDegree<3, 0, 1> {
  static constexpr Int value = 14;
};

template <>
struct AvgDegree<3, 0, 3> {
  static constexpr Int value = 24;
};

}  // end namespace Omega_h

#endif
