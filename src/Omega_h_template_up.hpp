#ifndef OMEGA_H_TEMPLATE_UP_HPP
#define OMEGA_H_TEMPLATE_UP_HPP

/*! \file Omega_h_template_up.hpp
  \brief Defines the TemplateUp class.
  */

#include <Omega_h_defines.hpp>

namespace Omega_h {

/*! \details Describes the upward relationship between an entity of dimension
   (d+1) (called the upward entity) adjacent to an entity of dimension d (called
   the downward entity), both of which are part of the local boundary of a third
   entity of dimension > d+1 (called the parent entity).
 */
struct TemplateUp {
  /*! \brief The parent-local index of the upward entity */
  Int up;
  /*! \brief The upward-local index of the downward entity */
  Int which_down;
  /*! \brief Whether the downward entity is flipped with respect to its
             canonical direction in the upward entity */
  bool is_flipped;
};

struct SplitVertex {
  Int dim;
  Int which_down;
};

}  // namespace Omega_h

#endif
