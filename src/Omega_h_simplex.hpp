#ifndef OMEGA_H_SIMPLEX_HPP
#define OMEGA_H_SIMPLEX_HPP

#include <Omega_h_template_up.hpp>
#include <Omega_h_kokkos.hpp>

/*! \file Omega_h_simplex.hpp
  \brief Describes the canonical local boundary connectivity
         orderings for a single simplex, and other related properties
  \details A simplex is a generalization of a triangle in arbitrary dimensions.
           In Omega_h the relevant simplices are vertices, edges, triangles, and tetrahedra.
  */

namespace Omega_h {

/* TODO: make these constexpr, either with C++14 or lots of
   ternary operators */

/*! \brief Relates bounding simplex vertices the parent simplex's vertices
  \param elem_dim The parent simplex's dimension
  \param bdry_dim The bounding simplex's dimension
  \param which_bdry The parent-local index of the bounding simplex
  \param which_vert The bounding-local index of the vertex
  \returns The parent-local index of the vertex
  */
OMEGA_H_INLINE Int simplex_down_template(
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
/* TODO: rename to singular_simplex_names */
extern char const* const singular_names[DIMS];
/* TODO: rename to plural_simplex_names */
extern char const* const plural_names[DIMS];

/* TODO: rename to SimplexAverageDegree */
template <Int dim, Int low, Int high>
struct AvgDegree;

template <>
struct AvgDegree<1, 0, 1> {
  static constexpr Int value = 2;
};

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

// this is some measure of the expected average quality
// of an element in a typical mesh.
// it can be derived by looking at meshes like the
// ones generated by osh_box, or just measuring average
// quality on real meshes.
// real measurements on the regression tests suggest
// the between 70% and 75% is reasonable for tetrahedra,
// although the full range goes down below 70% and above 80%
template <Int dim>
struct AvgQuality;

}  // end namespace Omega_h

#endif
