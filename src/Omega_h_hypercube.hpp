#ifndef OMEGA_H_HYPERCUBE_HPP
#define OMEGA_H_HYPERCUBE_HPP

#include <Omega_h_kokkos.hpp>
#include <Omega_h_template_up.hpp>

/*! \file Omega_h_hypercube.hpp
  \brief Describes the canonical local boundary connectivity
         orderings for a single hypercube, and other related properties
  \details A hypercube is a generalization of a square in arbitrary dimensions.
           In Omega_h the relevant hypercubes are vertices, edges,
  quadrilaterals, and hexahedra.
  */

namespace Omega_h {

/* TODO: make these constexpr, either with C++14 or lots of
   ternary operators */

/*! \brief Relates bounding hypercube vertices the parent hypercube's vertices
  \param elem_dim The parent hypercube's dimension
  \param bdry_dim The bounding hypercube's dimension
  \param which_bdry The parent-local index of the bounding hypercube
  \param which_vert The bounding-local index of the vertex
  \returns The parent-local index of the vertex
  */
OMEGA_H_INLINE Int hypercube_down_template(
    Int elem_dim, Int bdry_dim, Int which_bdry, Int which_vert) {
  switch (elem_dim) {
    case 1:
      switch (bdry_dim) {
        case 0:
          return which_bdry;
      }
      return -1;
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
	      return -1;
            case 1:
              switch (which_vert) {
                case 0:
                  return 1;
                case 1:
                  return 2;
              }
	      return -1;
            case 2:
              switch (which_vert) {
                case 0:
                  return 2;
                case 1:
                  return 3;
              }
	      return -1;
            case 3:
              switch (which_vert) {
                case 0:
                  return 3;
                case 1:
                  return 0;
              }
	      return -1;
          }
      }
      return -1;
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
	      return -1;
            case 1:
              switch (which_vert) {
                case 0:
                  return 1;
                case 1:
                  return 2;
              }
	      return -1;
            case 2:
              switch (which_vert) {
                case 0:
                  return 2;
                case 1:
                  return 3;
              }
	      return -1;
            case 3:
              switch (which_vert) {
                case 0:
                  return 3;
                case 1:
                  return 0;
              }
	      return -1;
            case 4:
              switch (which_vert) {
                case 0:
                  return 0;
                case 1:
                  return 4;
              }
	      return -1;
            case 5:
              switch (which_vert) {
                case 0:
                  return 1;
                case 1:
                  return 5;
              }
	      return -1;
            case 6:
              switch (which_vert) {
                case 0:
                  return 2;
                case 1:
                  return 6;
              }
	      return -1;
            case 7:
              switch (which_vert) {
                case 0:
                  return 3;
                case 1:
                  return 7;
              }
	      return -1;
            case 8:
              switch (which_vert) {
                case 0:
                  return 4;
                case 1:
                  return 5;
              }
	      return -1;
            case 9:
              switch (which_vert) {
                case 0:
                  return 5;
                case 1:
                  return 6;
              }
	      return -1;
            case 10:
              switch (which_vert) {
                case 0:
                  return 6;
                case 1:
                  return 7;
              }
	      return -1;
            case 11:
              switch (which_vert) {
                case 0:
                  return 7;
                case 1:
                  return 4;
              }
	      return -1;
          }
	  return -1;
        case 2:
          switch (which_bdry) {
            case 0:
              switch (which_vert) {
                case 0:
                  return 0;
                case 1:
                  return 3;
                case 2:
                  return 2;
                case 3:
                  return 1;
              }
	      return -1;
            case 1:
              switch (which_vert) {
                case 0:
                  return 0;
                case 1:
                  return 1;
                case 2:
                  return 5;
                case 3:
                  return 4;
              }
	      return -1;
            case 2:
              switch (which_vert) {
                case 0:
                  return 1;
                case 1:
                  return 2;
                case 2:
                  return 6;
                case 3:
                  return 5;
              }
	      return -1;
            case 3:
              switch (which_vert) {
                case 0:
                  return 2;
                case 1:
                  return 3;
                case 2:
                  return 7;
                case 3:
                  return 6;
              }
	      return -1;
            case 4:
              switch (which_vert) {
                case 0:
                  return 3;
                case 1:
                  return 0;
                case 2:
                  return 4;
                case 3:
                  return 7;
              }
	      return -1;
            case 5:
              switch (which_vert) {
                case 0:
                  return 4;
                case 1:
                  return 5;
                case 2:
                  return 6;
                case 3:
                  return 7;
              }
	      return -1;
          }
	  return -1;
      }
      return -1;
  }
  return -1;
}

OMEGA_H_INLINE TemplateUp hypercube_up_template(
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
                  return {3, 1, 0};
                case 2:
                  return {4, 0, 0};
              }
              return {-1, -1, 0};
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0};
                case 1:
                  return {0, 1, 0};
                case 2:
                  return {5, 0, 0};
              }
              return {-1, -1, 0};
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0};
                case 1:
                  return {1, 1, 0};
                case 2:
                  return {6, 0, 0};
              }
              return {-1, -1, 0};
            case 3:
              switch (which_up) {
                case 0:
                  return {3, 0, 0};
                case 1:
                  return {2, 1, 0};
                case 2:
                  return {7, 0, 0};
              }
              return {-1, -1, 0};
            case 4:
              switch (which_up) {
                case 0:
                  return {8, 0, 0};
                case 1:
                  return {11, 1, 0};
                case 2:
                  return {4, 1, 0};
              }
              return {-1, -1, 0};
            case 5:
              switch (which_up) {
                case 0:
                  return {9, 0, 0};
                case 1:
                  return {8, 1, 0};
                case 2:
                  return {5, 1, 0};
              }
              return {-1, -1, 0};
            case 6:
              switch (which_up) {
                case 0:
                  return {10, 0, 0};
                case 1:
                  return {9, 1, 0};
                case 2:
                  return {6, 1, 0};
              }
              return {-1, -1, 0};
            case 7:
              switch (which_up) {
                case 0:
                  return {11, 0, 0};
                case 1:
                  return {10, 1, 0};
                case 2:
                  return {7, 1, 0};
              }
              return {-1, -1, 0};
          }
	  return {-1, -1, 0};
        case 1:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 3, 1};
                case 1:
                  return {1, 0, 0};
              }
              return {-1, -1, 0};
            case 1:
              switch (which_up) {
                case 0:
                  return {0, 2, 1};
                case 1:
                  return {2, 0, 0};
              }
              return {-1, -1, 0};
            case 2:
              switch (which_up) {
                case 0:
                  return {0, 1, 1};
                case 1:
                  return {3, 0, 0};
              }
              return {-1, -1, 0};
            case 3:
              switch (which_up) {
                case 0:
                  return {0, 0, 1};
                case 1:
                  return {4, 0, 0};
              }
              return {-1, -1, 0};
            case 4:
              switch (which_up) {
                case 0:
                  return {1, 3, 1};
                case 1:
                  return {4, 1, 0};
              }
              return {-1, -1, 0};
            case 5:
              switch (which_up) {
                case 0:
                  return {2, 3, 1};
                case 1:
                  return {1, 1, 0};
              }
              return {-1, -1, 0};
            case 6:
              switch (which_up) {
                case 0:
                  return {3, 3, 1};
                case 1:
                  return {2, 1, 0};
              }
              return {-1, -1, 0};
            case 7:
              switch (which_up) {
                case 0:
                  return {4, 3, 1};
                case 1:
                  return {3, 1, 0};
              }
              return {-1, -1, 0};
            case 8:
              switch (which_up) {
                case 0:
                  return {5, 0, 0};
                case 1:
                  return {1, 2, 1};
              }
              return {-1, -1, 0};
            case 9:
              switch (which_up) {
                case 0:
                  return {5, 1, 0};
                case 1:
                  return {2, 2, 1};
              }
              return {-1, -1, 0};
            case 10:
              switch (which_up) {
                case 0:
                  return {5, 2, 0};
                case 1:
                  return {3, 2, 1};
              }
              return {-1, -1, 0};
            case 11:
              switch (which_up) {
                case 0:
                  return {5, 3, 0};
                case 1:
                  return {4, 2, 1};
              }
              return {-1, -1, 0};
          }
      }
      return {-1, -1, 0};
    case 2:
      switch (bdry_dim) {
        case 0:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {3, 1, 0};
                case 1:
                  return {0, 0, 0};
              }
              return {-1, -1, 0};
            case 1:
              switch (which_up) {
                case 0:
                  return {0, 1, 0};
                case 1:
                  return {1, 0, 0};
              }
              return {-1, -1, 0};
            case 2:
              switch (which_up) {
                case 0:
                  return {1, 1, 0};
                case 1:
                  return {2, 0, 0};
              }
              return {-1, -1, 0};
            case 3:
              switch (which_up) {
                case 0:
                  return {2, 1, 0};
                case 1:
                  return {3, 0, 0};
              }
              return {-1, -1, 0};
          }
      }
  }
  return {-1, -1, 0};
};

OMEGA_H_INLINE constexpr Int hypercube_degree(Int from_dim, Int to_dim) {
  // clang-format off
  return (from_dim == 0 ?
           (to_dim == 0 ? 1 : -1) :
         (from_dim == 1 ?
           (to_dim == 0 ? 2 :
           (to_dim == 1 ? 1 : -1)) :
         (from_dim == 2 ?
           (to_dim == 0 ? 4 :
           (to_dim == 1 ? 4 :
           (to_dim == 2 ? 1 : -1))) :
         (from_dim == 3 ?
           (to_dim == 0 ? 8 :
           (to_dim == 1 ? 12 :
           (to_dim == 2 ? 6 :
           (to_dim == 3 ? 1 : -1)))) : -1))));
  // clang-format on
}

constexpr char const* hypercube_singular_name(Int dim) {
  return (dim == 3 ? "hex"
                   : (dim == 2 ? "quad"
                               : (dim == 1 ? "edge"
                                           : (dim == 0 ? "vertex" : nullptr))));
}

constexpr char const* hypercube_plural_name(Int dim) {
  return (dim == 3
              ? "hexes"
              : (dim == 2 ? "quads"
                          : (dim == 1 ? "edges"
                                      : (dim == 0 ? "vertices" : nullptr))));
}

// every interior split entity can be seen as the dual of a corresponding
// boundary entity
OMEGA_H_INLINE constexpr Int hypercube_split_degree(
    Int parent_dim, Int child_dim) {
  return hypercube_degree(parent_dim, parent_dim - child_dim);
}

OMEGA_H_INLINE SplitVertex hypercube_split_template(
    Int parent_dim, Int child_dim, Int which_child, Int which_child_vtx) {
  switch (parent_dim) {
    case 1:
      switch (child_dim) {
        case 0:
          return {1, 0};
        case 1:
          switch (which_child) {
            case 0:
              switch (which_child_vtx) {
                case 0:
                  return {0, 0};
                case 1:
                  return {1, 0};
              }
              return {-1, -1};
            case 1:
              switch (which_child_vtx) {
                case 0:
                  return {1, 0};
                case 1:
                  return {0, 1};
              }
              return {-1, -1};
          }
      }
    case 2:
      switch (child_dim) {
        case 0:
          return {2, 0};
        case 1:
          switch (which_child_vtx) {
            case 0:
              return {1, which_child};
            case 1:
              return {2, 0};
          }
          return {-1, -1};
        case 2:
          switch (which_child) {
            case 0:
              switch (which_child_vtx) {
                case 0:
                  return {0, 0};
                case 1:
                  return {1, 0};
                case 2:
                  return {2, 0};
                case 3:
                  return {1, 3};
              }
              return {-1, -1};
            case 1:
              switch (which_child_vtx) {
                case 0:
                  return {1, 0};
                case 1:
                  return {0, 1};
                case 2:
                  return {1, 1};
                case 3:
                  return {2, 0};
              }
              return {-1, -1};
            case 2:
              switch (which_child_vtx) {
                case 0:
                  return {2, 0};
                case 1:
                  return {1, 1};
                case 2:
                  return {0, 2};
                case 3:
                  return {1, 2};
              }
              return {-1, -1};
            case 3:
              switch (which_child_vtx) {
                case 0:
                  return {1, 3};
                case 1:
                  return {2, 0};
                case 2:
                  return {1, 2};
                case 3:
                  return {0, 3};
              }
              return {-1, -1};
          }
      }
    case 3:
      switch (child_dim) {
        case 0:
          return {3, 0};
        case 1:
          switch (which_child_vtx) {
            case 0:
              return {2, which_child};
            case 1:
              return {3, 0};
          }
          return {-1, -1};
        case 2:
          switch (which_child) {
            case 0:
              switch (which_child_vtx) {
                case 0:
                  return {1, 0};
                case 1:
                  return {2, 0};
                case 2:
                  return {3, 0};
                case 3:
                  return {2, 1};
              }
              return {-1, -1};
            case 1:
              switch (which_child_vtx) {
                case 0:
                  return {2, 0};
                case 1:
                  return {1, 1};
                case 2:
                  return {2, 2};
                case 3:
                  return {3, 0};
              }
              return {-1, -1};
            case 2:
              switch (which_child_vtx) {
                case 0:
                  return {2, 0};
                case 1:
                  return {1, 2};
                case 2:
                  return {2, 3};
                case 3:
                  return {3, 0};
              }
              return {-1, -1};
            case 3:
              switch (which_child_vtx) {
                case 0:
                  return {1, 3};
                case 1:
                  return {2, 0};
                case 2:
                  return {3, 0};
                case 3:
                  return {2, 4};
              }
              return {-1, -1};
            case 4:
              switch (which_child_vtx) {
                case 0:
                  return {1, 4};
                case 1:
                  return {2, 1};
                case 2:
                  return {3, 0};
                case 3:
                  return {2, 4};
              }
              return {-1, -1};
            case 5:
              switch (which_child_vtx) {
                case 0:
                  return {2, 1};
                case 1:
                  return {1, 5};
                case 2:
                  return {2, 2};
                case 3:
                  return {3, 0};
              }
              return {-1, -1};
            case 6:
              switch (which_child_vtx) {
                case 0:
                  return {3, 0};
                case 1:
                  return {2, 2};
                case 2:
                  return {1, 6};
                case 3:
                  return {2, 3};
              }
              return {-1, -1};
            case 7:
              switch (which_child_vtx) {
                case 0:
                  return {2, 4};
                case 1:
                  return {3, 0};
                case 2:
                  return {2, 3};
                case 3:
                  return {1, 7};
              }
              return {-1, -1};
            case 8:
              switch (which_child_vtx) {
                case 0:
                  return {2, 1};
                case 1:
                  return {3, 0};
                case 2:
                  return {2, 5};
                case 3:
                  return {1, 8};
              }
              return {-1, -1};
            case 9:
              switch (which_child_vtx) {
                case 0:
                  return {3, 0};
                case 1:
                  return {2, 2};
                case 2:
                  return {1, 9};
                case 3:
                  return {2, 5};
              }
              return {-1, -1};
            case 10:
              switch (which_child_vtx) {
                case 0:
                  return {3, 0};
                case 1:
                  return {2, 3};
                case 2:
                  return {1, 10};
                case 3:
                  return {2, 5};
              }
              return {-1, -1};
            case 11:
              switch (which_child_vtx) {
                case 0:
                  return {2, 4};
                case 1:
                  return {3, 0};
                case 2:
                  return {2, 5};
                case 3:
                  return {1, 11};
              }
              return {-1, -1};
          }
        case 3:
          switch (which_child) {
            case 0:
              switch (which_child_vtx) {
                case 0:
                  return {0, 0};
                case 1:
                  return {1, 0};
                case 2:
                  return {2, 0};
                case 3:
                  return {1, 3};
                case 4:
                  return {1, 4};
                case 5:
                  return {2, 1};
                case 6:
                  return {3, 0};
                case 7:
                  return {2, 4};
              }
              return {-1, -1};
            case 1:
              switch (which_child_vtx) {
                case 0:
                  return {1, 0};
                case 1:
                  return {0, 1};
                case 2:
                  return {1, 1};
                case 3:
                  return {2, 0};
                case 4:
                  return {2, 1};
                case 5:
                  return {1, 5};
                case 6:
                  return {2, 2};
                case 7:
                  return {3, 0};
              }
              return {-1, -1};
            case 2:
              switch (which_child_vtx) {
                case 0:
                  return {2, 0};
                case 1:
                  return {1, 1};
                case 2:
                  return {0, 2};
                case 3:
                  return {1, 2};
                case 4:
                  return {3, 0};
                case 5:
                  return {2, 2};
                case 6:
                  return {1, 6};
                case 7:
                  return {2, 3};
              }
              return {-1, -1};
            case 3:
              switch (which_child_vtx) {
                case 0:
                  return {1, 3};
                case 1:
                  return {2, 0};
                case 2:
                  return {1, 2};
                case 3:
                  return {0, 3};
                case 4:
                  return {2, 4};
                case 5:
                  return {3, 0};
                case 6:
                  return {2, 3};
                case 7:
                  return {1, 7};
              }
              return {-1, -1};
            case 4:
              switch (which_child_vtx) {
                case 0:
                  return {1, 4};
                case 1:
                  return {2, 1};
                case 2:
                  return {3, 0};
                case 3:
                  return {2, 4};
                case 4:
                  return {0, 4};
                case 5:
                  return {1, 8};
                case 6:
                  return {2, 5};
                case 7:
                  return {1, 11};
              }
              return {-1, -1};
            case 5:
              switch (which_child_vtx) {
                case 0:
                  return {2, 1};
                case 1:
                  return {1, 5};
                case 2:
                  return {2, 2};
                case 3:
                  return {3, 0};
                case 4:
                  return {1, 8};
                case 5:
                  return {0, 5};
                case 6:
                  return {1, 9};
                case 7:
                  return {2, 5};
              }
              return {-1, -1};
            case 6:
              switch (which_child_vtx) {
                case 0:
                  return {3, 0};
                case 1:
                  return {2, 2};
                case 2:
                  return {1, 6};
                case 3:
                  return {2, 3};
                case 4:
                  return {2, 5};
                case 5:
                  return {1, 9};
                case 6:
                  return {0, 6};
                case 7:
                  return {1, 10};
              }
              return {-1, -1};
            case 7:
              switch (which_child_vtx) {
                case 0:
                  return {2, 4};
                case 1:
                  return {3, 0};
                case 2:
                  return {2, 3};
                case 3:
                  return {1, 7};
                case 4:
                  return {1, 11};
                case 5:
                  return {2, 5};
                case 6:
                  return {1, 10};
                case 7:
                  return {0, 7};
              }
              return {-1, -1};
          }
      }
  }
  return {-1, -1};
}

}  // end namespace Omega_h

#endif
