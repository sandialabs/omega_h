#ifndef MARK_HPP
#define MARK_HPP

#include "internal.hpp"

namespace osh {

Read<I8> mark_exposed_sides(Mesh* mesh);

Read<I8> mark_down(Mesh* mesh, Int high_dim, Int low_dim,
                   Read<I8> marked_highs);
Read<I8> mark_up(Mesh* mesh, Int low_dim, Int high_dim, Read<I8> low_marked);

Read<I8> mark_by_class_dim(Mesh* mesh, Int ent_dim, Int class_dim);
Read<I8> mark_by_class(Mesh* mesh, Int ent_dim, Int class_dim, I32 class_id);

Read<I8> mark_by_owner(Mesh* mesh, Int ent_dim, I32 rank);

Read<I8> mark_dual_layers(Mesh* mesh, Read<I8> marks, Int nlayers);

GO count_owned_marks(Mesh* mesh, Int ent_dim, Read<I8> marks);

Read<I8> mark_sliver_layers(Mesh* mesh, Real qual_ceil, Int nlayers);

} //end namespace osh

#endif
