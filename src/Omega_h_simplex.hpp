#ifndef OMEGA_H_SIMPLEX_HPP
#define OMEGA_H_SIMPLEX_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_kokkos.hpp>

namespace Omega_h {

/* TODO: make these constexpr, either with C++14 or lots of
   ternary operators */
OMEGA_H_INLINE Int down_template(
    Int elem_dim, Int bdry_dim, Int which_bdry, Int which_vert) {
  switch (elem_dim) {
    case 1:
      switch (bdry_dim) {
        case 0:
          return which_bdry;
      }
    case 2:
      switch (bdry_dim) {
        case 0:
          return which_bdry;
        case 1:
          switch (which_bdry) {
            case 0:
              switch (which_vert) {
                case 0:
                  return 0;
                case 1:
                  return 1;
              }
            case 1:
              switch (which_vert) {
                case 0:
                  return 1;
                case 1:
                  return 2;
              }
            case 2:
              switch (which_vert) {
                case 0:
                  return 2;
                case 1:
                  return 0;
              }
          }
      }
    case 3:
      switch (bdry_dim) {
        case 0:
          return which_bdry;
        case 1:
          switch (which_bdry) {
            case 0:
              switch (which_vert) {
                case 0:
                  return 0;
                case 1:
                  return 1;
              }
            case 1:
              switch (which_vert) {
                case 0:
                  return 1;
                case 1:
                  return 2;
              }
            case 2:
              switch (which_vert) {
                case 0:
                  return 2;
                case 1:
                  return 0;
              }
            case 3:
              switch (which_vert) {
                case 0:
                  return 0;
                case 1:
                  return 3;
              }
            case 4:
              switch (which_vert) {
                case 0:
                  return 1;
                case 1:
                  return 3;
              }
            case 5:
              switch (which_vert) {
                case 0:
                  return 2;
                case 1:
                  return 3;
              }
          }
        case 2:
          switch (which_bdry) {
            case 0:
              switch (which_vert) {
                case 0:
                  return 0;
                case 1:
                  return 2;
                case 2:
                  return 1;
              }
            case 1:
              switch (which_vert) {
                case 0:
                  return 0;
                case 1:
                  return 1;
                case 2:
                  return 3;
              }
            case 2:
              switch (which_vert) {
                case 0:
                  return 1;
                case 1:
                  return 2;
                case 2:
                  return 3;
              }
            case 3:
              switch (which_vert) {
                case 0:
                  return 2;
                case 1:
                  return 0;
                case 2:
                  return 3;
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

OMEGA_H_INLINE TemplateUp up_template(
    Int elem_dim, Int bdry_dim, Int which_bdry, Int which_up) {
  switch (elem_dim) {
    case 3:
      switch (bdry_dim) {
        case 0:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 0, 0};
                case 1:
                  return {2, 1, 0};
                case 2:
                  return {3, 0, 0};
              }
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0};
                case 1:
                  return {0, 1, 0};
                case 2:
                  return {4, 0, 0};
              }
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0};
                case 1:
                  return {1, 1, 0};
                case 2:
                  return {5, 0, 0};
              }
            case 3:
              switch (which_up) {
                case 0:
                  return {5, 1, 0};
                case 1:
                  return {4, 1, 0};
                case 2:
                  return {3, 1, 0};
              }
          }
        case 1:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 2, 1};
                case 1:
                  return {1, 0, 0};
              }
            case 1:
              switch (which_up) {
                case 0:
                  return {0, 1, 1};
                case 1:
                  return {2, 0, 0};
              }
            case 2:
              switch (which_up) {
                case 0:
                  return {0, 0, 1};
                case 1:
                  return {3, 0, 0};
              }
            case 3:
              switch (which_up) {
                case 0:
                  return {1, 2, 1};
                case 1:
                  return {3, 1, 0};
              }
            case 4:
              switch (which_up) {
                case 0:
                  return {2, 2, 1};
                case 1:
                  return {1, 1, 0};
              }
            case 5:
              switch (which_up) {
                case 0:
                  return {3, 2, 1};
                case 1:
                  return {2, 1, 0};
              }
          }
      }
    case 2:
      switch (bdry_dim) {
        case 0:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 0, 0};
                case 1:
                  return {2, 1, 0};
              }
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0};
                case 1:
                  return {0, 1, 0};
              }
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0};
                case 1:
                  return {1, 1, 0};
              }
          }
      }
  }
  return {-1, -1, true};
};

OMEGA_H_INLINE Int opposite_template(Int elem_dim, Int bdry_dim, Int which_bdry) {
  switch (elem_dim) {
    case 3:
      switch (bdry_dim) {
        case 0:
          switch (which_bdry) {
            case 0:
              return 2;
            case 1:
              return 3;
            case 2:
              return 1;
            case 3:
              return 0;
          }
        case 1:
          switch (which_bdry) {
            case 0:
              return 5;
            case 1:
              return 3;
            case 2:
              return 4;
            case 3:
              return 1;
            case 4:
              return 2;
            case 5:
              return 0;
          }
        case 2:
          switch (which_bdry) {
            case 0:
              return 3;
            case 1:
              return 2;
            case 2:
              return 0;
            case 3:
              return 1;
          }
      }
    case 2:
      switch (bdry_dim) {
        case 0:
          switch (which_bdry) {
            case 0:
              return 1;
            case 1:
              return 2;
            case 2:
              return 0;
          }
        case 1:
          switch (which_bdry) {
            case 0:
              return 2;
            case 1:
              return 0;
            case 2:
              return 1;
          }
      }
    case 1:
      return (1 - which_bdry);
  }
  return -1;
}

OMEGA_H_INLINE Int simplex_degree(Int from_dim, Int to_dim) {
  switch (from_dim) {
    case 0:
      return 1;
    case 1:
      switch (to_dim) {
        case 0:
          return 2;
        case 1:
          return 1;
      }
    case 2:
      switch (to_dim) {
        case 0:
          return 3;
        case 1:
          return 3;
        case 2:
          return 1;
      }
    case 3:
      switch (to_dim) {
        case 0:
          return 4;
        case 1:
          return 6;
        case 2:
          return 4;
        case 3:
          return 1;
      }
  }
  return -1;
}

/* TODO: replace simplex_degrees with simplex_degree() */
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
