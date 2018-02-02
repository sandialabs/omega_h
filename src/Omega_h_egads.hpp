#ifndef OMEGA_H_EGADS_HPP
#define OMEGA_H_EGADS_HPP

#include <Omega_h_array.hpp>
#include <string>

namespace Omega_h {

class Mesh;
struct Egads;

Egads* egads_load(std::string const& filename);
void egads_classify(Egads* eg, int nadj_faces, int const adj_face_ids[],
    int* class_dim, int* class_id);
void egads_free(Egads* eg);

void egads_reclassify(Mesh* mesh, Egads* eg);
Reals egads_get_snap_warp(Mesh* mesh, Egads* eg, bool verbose);

}  // namespace Omega_h

#endif
