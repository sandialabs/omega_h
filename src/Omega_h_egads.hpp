#ifndef OMEGA_H_EGADS_HPP
#define OMEGA_H_EGADS_HPP

#include <string>

namespace Omega_h {

struct Egads;

Egads* egads_load(std::string const& filename);
void egads_classify(Egads* eg, int nadj_faces, int const adj_face_ids[],
    int* class_dim, int* class_id);
void egads_free(Egads* eg);

}

#endif
