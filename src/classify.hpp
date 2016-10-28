#ifndef CLASSIFY_HPP
#define CLASSIFY_HPP

#include "internal.hpp"

namespace Omega_h {

void classify_sides_by_exposure(Mesh* mesh, Read<I8> side_is_exposed);
void classify_hinges_by_sharpness(
    Mesh* mesh, Read<I8> hinge_is_exposed, Read<I8> hinge_is_sharp);
void classify_elements(Mesh* mesh);

/* this function is meant to take in any amount
 * of existing classification and do its best
 * to derive as much of the classification for
 * the rest of the mesh as possible.
 */
void finalize_classification(Mesh* mesh);

/* given a set of equal-order entities
 * (entities with the same dimension as the
 *  model entity they are classified on)
 * defined by their vertices, this function
 * will set the classification dimensions and IDs
 * for all entities of that dimension.
 * this function is typically called prior
 * to using finalize_classification()
 */
void classify_equal_order(Mesh* mesh, Int ent_dim, LOs eqv2v,
    Read<LO> eq_class_ids);

}  // end namespace Omega_h

#endif
