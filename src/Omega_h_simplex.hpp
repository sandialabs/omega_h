#ifndef OMEGA_H_SIMPLEX_HPP
#define OMEGA_H_SIMPLEX_HPP

#include <Omega_h_template_up.hpp>

/* Omega_h_simplex.hpp
  \brief Describes the canonical local boundary connectivity
         orderings for a single simplex, and other related properties
  \details A simplex is a generalization of a triangle in arbitrary dimensions.
           In Omega_h the relevant simplices are vertices, edges, triangles, and
  tetrahedra.
  */

namespace Omega_h {

/*! \brief Relates bounding simplex vertices to the parent simplex's vertices
  \param elem_dim The parent simplex's dimension
  \param bdry_dim The bounding simplex's dimension
  \param which_bdry The parent-local index of the bounding simplex
  \param which_vert The bounding-local index of the vertex
  \returns The parent-local index of the vertex
  */
constexpr OMEGA_H_INLINE Int simplex_down_template(
    Int elem_dim, Int bdry_dim, Int which_bdry, Int which_vert) {
  // clang-format off
  return (elem_dim == 3 ?
           (bdry_dim == 2 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 2 :
               (which_vert == 2 ? 1 : -1))) :
             (which_bdry == 1 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 1 :
               (which_vert == 2 ? 3 : -1))) :
             (which_bdry == 2 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 2 :
               (which_vert == 2 ? 3 : -1))) :
             (which_bdry == 3 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 0 :
               (which_vert == 2 ? 3 : -1))) : -1)))) :
           (bdry_dim == 1 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 1 : -1)) :
             (which_bdry == 1 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 2 : -1)) :
             (which_bdry == 2 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 0 : -1)) :
             (which_bdry == 3 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 3 : -1)) :
             (which_bdry == 4 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 3 : -1)) :
             (which_bdry == 5 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 3 : -1)) : -1)))))) :
           (bdry_dim == 0 ?
             (which_bdry == 0 ? 0 :
             (which_bdry == 1 ? 1 :
             (which_bdry == 2 ? 2 :
             (which_bdry == 3 ? 3 : -1)))) : -1))) :
         (elem_dim == 2 ?
           (bdry_dim == 1 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 1 : -1)) :
             (which_bdry == 1 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 2 : -1)) :
             (which_bdry == 2 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 0 : -1)) : -1))) :
           (bdry_dim == 0 ?
             (which_bdry == 0 ? 0 :
             (which_bdry == 1 ? 1 :
             (which_bdry == 2 ? 2 : -1))) : -1)) :
         (elem_dim == 1 ?
           (bdry_dim == 0 ?
             (which_bdry == 0 ? 0 :
             (which_bdry == 1 ? 1 : -1)) : -1) : -1)));
  // clang-format on
}

/* TODO: make these constexpr, either with C++14 or lots of
   ternary operators */

OMEGA_H_INLINE TemplateUp simplex_up_template(
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
              return {-1, -1, 0};
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0};
                case 1:
                  return {0, 1, 0};
                case 2:
                  return {4, 0, 0};
              }
              return {-1, -1, 0};
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0};
                case 1:
                  return {1, 1, 0};
                case 2:
                  return {5, 0, 0};
              }
              return {-1, -1, 0};
            case 3:
              switch (which_up) {
                case 0:
                  return {5, 1, 0};
                case 1:
                  return {4, 1, 0};
                case 2:
                  return {3, 1, 0};
              }
              return {-1, -1, 0};
          }
          return {-1, -1, 0};
        case 1:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 2, 1};
                case 1:
                  return {1, 0, 0};
              }
              return {-1, -1, 0};
            case 1:
              switch (which_up) {
                case 0:
                  return {0, 1, 1};
                case 1:
                  return {2, 0, 0};
              }
              return {-1, -1, 0};
            case 2:
              switch (which_up) {
                case 0:
                  return {0, 0, 1};
                case 1:
                  return {3, 0, 0};
              }
              return {-1, -1, 0};
            case 3:
              switch (which_up) {
                case 0:
                  return {1, 2, 1};
                case 1:
                  return {3, 1, 0};
              }
              return {-1, -1, 0};
            case 4:
              switch (which_up) {
                case 0:
                  return {2, 2, 1};
                case 1:
                  return {1, 1, 0};
              }
              return {-1, -1, 0};
            case 5:
              switch (which_up) {
                case 0:
                  return {3, 2, 1};
                case 1:
                  return {2, 1, 0};
              }
              return {-1, -1, 0};
          }
          return {-1, -1, 0};
      }
      return {-1, -1, 0};
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
              return {-1, -1, 0};
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0};
                case 1:
                  return {0, 1, 0};
              }
              return {-1, -1, 0};
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0};
                case 1:
                  return {1, 1, 0};
              }
              return {-1, -1, 0};
          }
          return {-1, -1, 0};
      }
      return {-1, -1, 0};
  }
  return {-1, -1, 0};
};

OMEGA_H_INLINE Int simplex_opposite_template(
    Int elem_dim, Int bdry_dim, Int which_bdry) {
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
          return -1;
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
          return -1;
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
          return -1;
      }
      return -1;
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
          return -1;
        case 1:
          switch (which_bdry) {
            case 0:
              return 2;
            case 1:
              return 0;
            case 2:
              return 1;
          }
          return -1;
      }
      return -1;
    case 1:
      switch (bdry_dim) {
        case 0:
          return (1 - which_bdry);
      }
      return -1;
  }
  return -1;
}

OMEGA_H_INLINE constexpr Int simplex_degree(Int from_dim, Int to_dim) {
  // clang-format off
  return (from_dim == 0 ? 1 :
         (from_dim == 1 ?
           (to_dim == 0 ? 2 :
           (to_dim == 1 ? 1 : -1)) :
         (from_dim == 2 ?
           (to_dim == 0 ? 3 :
           (to_dim == 1 ? 3 :
           (to_dim == 2 ? 1 : -1))) :
         (from_dim == 3 ?
           (to_dim == 0 ? 4 :
           (to_dim == 1 ? 6 :
           (to_dim == 2 ? 4 :
           (to_dim == 3 ? 1 : -1)))) : -1))));
  // clang-format on
}

constexpr char const* simplex_singular_name(Int dim) {
  return (dim == 3 ? "tet"
                   : (dim == 2 ? "triangle"
                               : (dim == 1 ? "edge"
                                           : (dim == 0 ? "vertex" : nullptr))));
}

constexpr char const* simplex_plural_name(Int dim) {
  return (dim == 3
              ? "tets"
              : (dim == 2 ? "triangles"
                          : (dim == 1 ? "edges"
                                      : (dim == 0 ? "vertices" : nullptr))));
}

template <Int dim, Int low, Int high>
struct SimplexAvgDegree;

template <>
struct SimplexAvgDegree<1, 0, 1> {
  static constexpr Int value = 2;
};

template <>
struct SimplexAvgDegree<2, 0, 1> {
  static constexpr Int value = 6;
};

template <>
struct SimplexAvgDegree<2, 0, 2> {
  static constexpr Int value = 6;
};

template <>
struct SimplexAvgDegree<3, 0, 1> {
  static constexpr Int value = 14;
};

template <>
struct SimplexAvgDegree<3, 0, 3> {
  static constexpr Int value = 24;
};

}  // end namespace Omega_h

#endif
