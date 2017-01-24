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

INLINE TemplateUp up_template(Int elem_dim, Int bdry_dim, Int which_bdry, Int which_up) {
  switch (elem_dim) {
    case 3:
      switch (bdry_dim) {
        case 0:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0: return {0,0,0};
                case 1: return {2,1,0};
                case 2: return {3,0,0};
              }
            case 1:
              switch (which_up) {
                case 0: return {1,0,0};
                case 1: return {0,1,0};
                case 2: return {4,0,0};
              }
            case 2:
              switch (which_up) {
                case 0: return {2,0,0};
                case 1: return {1,1,0};
                case 2: return {5,0,0};
              }
            case 3:
              switch (which_up) {
                case 0: return {5,1,0};
                case 1: return {4,1,0};
                case 2: return {3,1,0};
              }
          }
        case 1:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0: return {0,2,1};
                case 1: return {1,0,0};
              }
            case 1:
              switch (which_up) {
                case 0: return {0,1,1};
                case 1: return {2,0,0};
              }
            case 2:
              switch (which_up) {
                case 0: return {0,0,1};
                case 1: return {3,0,0};
              }
            case 3:
              switch (which_up) {
                case 0: return {1,2,1};
                case 1: return {3,1,0};
              }
            case 4:
              switch (which_up) {
                case 0: return {2,2,1};
                case 1: return {1,1,0};
              }
            case 5:
              switch (which_up) {
                case 0: return {3,2,1};
                case 1: return {2,1,0};
              }
          }
      }
    case 2:
      switch (bdry_dim) {
        case 0:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0: return {0,0,0};
                case 1: return {2,1,0};
              }
            case 1:
              switch (which_up) {
                case 0: return {1,0,0};
                case 1: return {0,1,0};
              }
            case 2:
              switch (which_up) {
                case 0: return {2,0,0};
                case 1: return {1,1,0};
              }
          }
      }
  }
  return {-1,-1,-1};
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
