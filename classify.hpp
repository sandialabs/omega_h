#ifndef CLASSIFY_HPP
#define CLASSIFY_HPP

#include "internal.hpp"

namespace osh {

void classify_sides_by_exposure(Mesh* mesh, Read<I8> side_is_exposed);
void classify_hinges_by_sharpness(Mesh* mesh, Read<I8> hinge_is_exposed,
                                  Read<I8> hinge_is_sharp);
void classify_vertices_by_sharp_edges(Mesh* mesh, Read<I8> vert_is_exposed,
                                      Read<I8> edge_is_sharp);
void classify_elements(Mesh* mesh);
void classify_by_angles(Mesh* mesh, Real sharp_angle);
void project_classification(Mesh* mesh);

} //end namespace osh

#endif
