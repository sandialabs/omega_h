#include "Omega_h_file.hpp"

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"

#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"

namespace Omega_h {

namespace meshsim {

namespace {

void read_internal(pParMesh sm, Mesh* mesh) {
  (void)mesh;
  auto i = 0;
  pMesh m = PM_mesh(sm, 0);
  EIter edges = M_edgeIter(m);
  pEntity ent;
  while ((ent = EIter_next(edges))) {
    i++;
  }
  EIter_delete(edges);
  fprintf(stderr, "edges %d\n", i);
}

}  // end anonymous namespace

Mesh read(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm) {
  SimPartitionedMesh_start(NULL,NULL);
  SimModel_start();
  Sim_readLicenseFile(NULL);
  pNativeModel nm = NULL;
  pProgress p = NULL;
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  pParMesh sm = PM_load(mesh_fname.c_str(), g, p);
  auto mesh = Mesh(comm->library());
  meshsim::read_internal(sm, &mesh);
  M_release(sm);
  GM_release(g);
  SimModel_stop();
  SimPartitionedMesh_stop();
  return mesh;
}

}  // namespace meshsim

}  // end namespace Omega_h
