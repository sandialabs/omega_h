#ifndef INDSET_HPP
#define INDSET_HPP

namespace osh {

Read<I8> find_indset(Mesh* mesh, Int ent_dim, Reals quality,
                     Read<I8> candidates);

} //end namespace osh

#endif
