#include "consistent.hpp"

#include <iostream>

#include "array.hpp"
#include "simplices.hpp"
#include "tag.hpp"

namespace Omega_h {

template <typename T>
static bool is_consistent(Mesh* mesh, Int dim, Read<T> copy_data, Int ncomps) {
  auto synced_data = mesh->sync_array(dim, copy_data, ncomps);
  auto local_ok = (copy_data == synced_data);
  auto global_ok = mesh->comm()->reduce_and(local_ok);
  return global_ok;
}

template <typename T>
static bool is_consistent(Mesh* mesh, Int dim, TagBase const* tagbase) {
  auto tag = to<T>(tagbase);
  return is_consistent(mesh, dim, tag->array(), tagbase->ncomps());
}

static bool tags_are_consistent(Mesh* mesh, Int dim) {
  for (Int i = 0; i < mesh->ntags(dim); ++i) {
    auto tagbase = mesh->get_tag(dim, i);
    bool ok = false;
    switch (tagbase->type()) {
      case OMEGA_H_I8:
        ok = is_consistent<I8>(mesh, dim, tagbase);
        break;
      case OMEGA_H_I32:
        ok = is_consistent<I32>(mesh, dim, tagbase);
        break;
      case OMEGA_H_I64:
        ok = is_consistent<I64>(mesh, dim, tagbase);
        break;
      case OMEGA_H_F64:
        ok = is_consistent<Real>(mesh, dim, tagbase);
        break;
    }
    if (!ok && mesh->comm()->rank() == 0) {
      std::cerr << singular_names[dim] << " tag " << tagbase->name()
                << " is not consistent\n";
      return false;
    }
  }
  return true;
}

bool tags_are_consistent(Mesh* mesh) {
  for (Int dim = 0; dim <= mesh->dim(); ++dim) {
    if (!tags_are_consistent(mesh, dim)) return false;
  }
  return true;
}

}  // end namespace Omega_h
