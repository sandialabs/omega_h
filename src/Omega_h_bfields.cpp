#include "Omega_h_bfields.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_ghost.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_atomics.hpp"

namespace Omega_h {

#ifdef OMEGA_H_USE_CUDA
__host__
#endif
    void
    assign(Mesh& a, Mesh const& b) {
  a = b;
}

#define OMEGA_H_INST(T)                                                        \
  template void Mesh::change_tagToBoundary<T>(                                 \
      Int ent_dim, Int ncomps, std::string const& name);                       \
  template void Mesh::change_tagToMesh<T>(                                     \
      Int ent_dim, Int ncomps, std::string const& name);                       \
  template Read<T> Mesh::get_boundaryField_array<T>(                           \
      Int dim, std::string const& name) const;                                 \
  template void Mesh::add_boundaryField<T>(                                    \
      Int dim, std::string const& name, Int ncomps);                           \
  template void Mesh::add_boundaryField<T>(                                    \
      Int dim, std::string const& name, Int ncomps, Read<T> array,             \
      bool internal);                                                          \
  template void Mesh::set_boundaryField_array(                                 \
      Int dim, std::string const& name, Read<T> array, bool internal);         
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

bool Mesh::has_revClass (Int edim) const {
  check_dim2 (edim);
  return bool (revClass_[edim]);
}

Adj Mesh::get_revClass (Int edim) const {
  check_dim2 (edim);
  OMEGA_H_CHECK (has_revClass (edim));
  return *(revClass_[edim]);
}

Adj Mesh::derive_revClass (Int edim) {
  OMEGA_H_TIME_FUNCTION;
  check_dim2 (edim);

  auto class_ids_all = get_array<ClassId>(edim, "class_id");
  auto class_dim = get_array<I8>(edim, "class_dim");

  Write<LO> class_ids_w(nents(edim), -1, "edim_classIds");
  auto edim_classid = OMEGA_H_LAMBDA (LO i) {
    if (class_dim[i] == edim) {
      class_ids_w[i] = class_ids_all[i];
    }
  };
  parallel_for(nents(edim), std::move(edim_classid));
  auto class_ids = LOs(class_ids_w);

  auto const n_gents = get_max(class_ids) + 1;
  Write<LO> degree(n_gents, 0, "rc_degrees");
  auto count_degree = OMEGA_H_LAMBDA (LO i) {
    if (class_ids[i] >= 0) {
      auto const gent_id = class_ids[i];
      atomic_increment(&degree[gent_id]);
    }
  };
  parallel_for(nents(edim), std::move(count_degree));
  auto a2ab_r = offset_scan(Read<LO>(degree), "rc_a2ab");

  auto const total_ments = get_sum(Read<LO>(degree));
  Write<LO> ab2b(total_ments, 0, "rc_ab2b");
  Write<LO> positions(n_gents, 0, "rc_positions");

  auto get_values = OMEGA_H_LAMBDA (LO i) {
    if (class_ids[i] >= 0) {
      auto const gent_id = class_ids[i];
      auto const first = a2ab_r[gent_id];
      auto const j = atomic_fetch_add(&positions[gent_id], 1);
      ab2b[first + j] = i;
    }
  };
  parallel_for(nents(edim), std::move(get_values));

  return Adj(a2ab_r, LOs(ab2b));

}

Adj Mesh::ask_revClass (Int edim) {
  OMEGA_H_TIME_FUNCTION;
  check_dim2 (edim);
  if (has_revClass (edim)) {
    return get_revClass (edim);
  }
  Adj derived_rc = derive_revClass (edim);
  revClass_[edim] = std::make_shared<Adj>(derived_rc);
  return derived_rc;
}

Adj Mesh::ask_revClass (Int edim, LOs g_ids) {
  OMEGA_H_TIME_FUNCTION;
  check_dim2 (edim);
  if (!g_ids.size()) {
    fprintf(stderr, "Model entity IDs cannot be empty\n");
    OMEGA_H_NORETURN(Adj());
  }
  auto edim_rc = ask_revClass(edim);
  auto n_gents = g_ids.size();
  auto max_gent_id = get_max(g_ids);
  auto max_gent_id_avail = (edim_rc.a2ab.size() - 1);
  if (max_gent_id > max_gent_id_avail) {
    fprintf(stderr, "input model id %d greater than local max. model id\n",
            max_gent_id);
    OMEGA_H_NORETURN(Adj());
  }

  auto ab2b = edim_rc.ab2b;
  auto a2ab = edim_rc.a2ab;
  Write<LO> degree(max_gent_id+1, 0, "new_rc_degrees");

  auto count = OMEGA_H_LAMBDA (LO i) {
    auto gent = g_ids[i];
    auto start = a2ab[gent];
    auto end = a2ab[gent + 1];
    degree[gent] = end - start;
  };
  parallel_for(n_gents, std::move(count));

  auto total_ments = get_sum(Read<LO>(degree));
  auto new_a2ab_r = offset_scan(Read<LO>(degree), "new_rc_a2ab");
  Write<LO> new_ab2b(total_ments, 0, "new_rc_ab2b");

  auto get_values = OMEGA_H_LAMBDA (LO i) {
    auto gent = g_ids[i];
    auto start = a2ab[gent];
    auto end = a2ab[gent + 1];
    for (LO j = start; j < end; ++j) {
      new_ab2b[new_a2ab_r[gent] + j - start] = ab2b[j];
    }
  };
  parallel_for(n_gents, std::move(get_values));

  return Adj(new_a2ab_r, LOs(new_ab2b));
}

template <typename T>
void Mesh::change_tagToBoundary(Int ent_dim, Int ncomps,
                                std::string const &name) {

  OMEGA_H_TIME_FUNCTION;
  auto mesh_field = get_array<T>(ent_dim, name);
  OMEGA_H_CHECK (mesh_field.size() == nents(ent_dim)*ncomps);

  auto boundary_ids = (ask_revClass(ent_dim)).ab2b;
  auto n_bEnts = boundary_ids.size();
  Write<T> b_field(n_bEnts*ncomps);
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto id = boundary_ids[i];
    for (LO n = 0; n < ncomps; ++n) {
      b_field[i*ncomps + n] = mesh_field[id*ncomps + n];
    }
  };
  parallel_for(n_bEnts, f, "get_bdryField");

  set_tag<T>(ent_dim, name, Read<T>(b_field));

  return;
}

template <typename T>
void Mesh::change_tagToMesh(Int ent_dim, Int ncomps,
                            std::string const &name) {

  OMEGA_H_TIME_FUNCTION;
  auto boundary_field = get_array<T>(ent_dim, name);
  auto n_ents = nents (ent_dim);
  OMEGA_H_CHECK (boundary_field.size() <= n_ents*ncomps);

  auto boundary_ids = (ask_revClass(ent_dim)).ab2b;
  auto n_bEnts = boundary_ids.size();

  Write<T> mesh_field (n_ents*ncomps, OMEGA_H_INTERIOR_VAL);

  auto f = OMEGA_H_LAMBDA (LO i) {
    auto id = boundary_ids[i];
    for (LO n = 0; n < ncomps; ++n) {
      if (mesh_field[id*ncomps + n] == OMEGA_H_INTERIOR_VAL) {
        mesh_field[id*ncomps + n] = boundary_field[i*ncomps + n];
      }
    }
  };
  parallel_for(n_bEnts, f, "get_fieldFromBdry");

  set_tag<T>(ent_dim, name, Read<T>(mesh_field));

  return;
}

template <typename T>
Read<T> Mesh::get_boundaryField_array
  (Int ent_dim, std::string const& name) const {

  size_t found = name.find("_boundary");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_boundary' at end of field name\n");
    OMEGA_H_NORETURN(Read<T>());
  }

  std::string new_name = name;
  new_name.append("_boundary");
  OMEGA_H_CHECK(has_tag(ent_dim, new_name));
  return get_tag<T>(ent_dim, new_name)->array();

}

template <typename T>
void Mesh::add_boundaryField(Int ent_dim, std::string const& name,
                             Int ncomps) {

  size_t found = name.find("_boundary");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_boundary' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_boundary");
  OMEGA_H_CHECK(!has_tag(ent_dim, new_name));
  add_tag<T>(ent_dim, new_name, ncomps);

  return;
}

template <typename T>
void Mesh::add_boundaryField(Int ent_dim, std::string const& name, Int ncomps,
  Read<T> array, bool internal) {

  size_t found = name.find("_boundary");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_boundary' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_boundary");
  OMEGA_H_CHECK(!has_tag(ent_dim, new_name));
  add_tag<T>(ent_dim, new_name, ncomps, array, internal);

  return;
}

bool Mesh::has_boundaryField(Int ent_dim, std::string const& name) const {

  size_t found = name.find("_boundary");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_boundary' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_boundary");
  if (!has_ents(ent_dim)) return false;
  return has_tag(ent_dim, new_name);

}

void Mesh::remove_boundaryField(Int ent_dim, std::string const& name) {

  size_t found = name.find("_boundary");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_boundary' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_boundary");
  if (!has_ents(ent_dim)) return;
  remove_tag(ent_dim, new_name);

}

template <typename T>
void Mesh::set_boundaryField_array
  (Int ent_dim, std::string const& name, Read<T> array, bool internal) {

  size_t found = name.find("_boundary");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_boundary' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_boundary");
  OMEGA_H_CHECK(has_tag(ent_dim, new_name));
  set_tag<T>(ent_dim, new_name, array, internal);

  return;
}

void Mesh::reduce_boundaryField(Int ent_dim, std::string const& name,
     Omega_h_Op op) {

  std::string new_name = name;
  new_name.append("_boundary");

  size_t found = name.find("_boundary");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_boundary' at end of field name\n");
  }
  else {
    OMEGA_H_CHECK(has_tag(ent_dim, new_name));
  }

  auto tagbase = get_tagbase(ent_dim, new_name);
  switch (tagbase->type()) {
    case OMEGA_H_I8: {

      if (nents(ent_dim))
        change_tagToMesh<I8> (ent_dim, tagbase->ncomps(), new_name);

      auto out = reduce_array(
          ent_dim, as<I8>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(ent_dim, new_name, out);

      change_tagToBoundary<I8> (ent_dim, tagbase->ncomps(), new_name);

      break;
    }
    case OMEGA_H_I32: {

      if (nents(ent_dim))
        change_tagToMesh<I32> (ent_dim, tagbase->ncomps(), new_name);

      auto out = reduce_array(
          ent_dim, as<I32>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(ent_dim, new_name, out);

      change_tagToBoundary<I32> (ent_dim, tagbase->ncomps(), new_name);

      break;
    }
    case OMEGA_H_I64: {

      if (nents(ent_dim)) 
        change_tagToMesh<I64> (ent_dim, tagbase->ncomps(), new_name);

      auto out = reduce_array(
          ent_dim, as<I64>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(ent_dim, new_name, out);

      change_tagToBoundary<I64> (ent_dim, tagbase->ncomps(), new_name);

      break;
    }
    case OMEGA_H_F64: {

      if (nents(ent_dim)) 
        change_tagToMesh<Real> (ent_dim, tagbase->ncomps(), new_name);

      auto out = reduce_array(
          ent_dim, as<Real>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(ent_dim, new_name, out);

      change_tagToBoundary<Real> (ent_dim, tagbase->ncomps(), new_name);

      break;
    }
  }
}

void Mesh::sync_boundaryField(Int ent_dim, std::string const& name) {

  std::string new_name = name;
  new_name.append("_boundary");

  size_t found = name.find("_boundary");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_boundary' at end of field name\n");
  }
  else {
    OMEGA_H_CHECK(has_tag(ent_dim, new_name));
  }

  auto tagbase = get_tagbase(ent_dim, new_name);
  switch (tagbase->type()) {
    case OMEGA_H_I8: {

      if (nents(ent_dim)) 
        change_tagToMesh<I8> (ent_dim, tagbase->ncomps(), new_name);

      auto out =
          sync_array(ent_dim, as<I8>(tagbase)->array(), tagbase->ncomps());
      set_tag(ent_dim, new_name, out);

      change_tagToBoundary<I8> (ent_dim, tagbase->ncomps(), new_name);

      break;
    }
    case OMEGA_H_I32: {

      if (nents(ent_dim)) 
        change_tagToMesh<I32> (ent_dim, tagbase->ncomps(), new_name);

      auto out =
          sync_array(ent_dim, as<I32>(tagbase)->array(), tagbase->ncomps());
      set_tag(ent_dim, new_name, out);

      change_tagToBoundary<I32> (ent_dim, tagbase->ncomps(), new_name);

      break;
    }
    case OMEGA_H_I64: {

      if (nents(ent_dim)) 
        change_tagToMesh<I64> (ent_dim, tagbase->ncomps(), new_name);

      auto out =
          sync_array(ent_dim, as<I64>(tagbase)->array(), tagbase->ncomps());
      set_tag(ent_dim, new_name, out);

      change_tagToBoundary<I64> (ent_dim, tagbase->ncomps(), new_name);

      break;
    }
    case OMEGA_H_F64: {

      if (nents(ent_dim))
        change_tagToMesh<Real> (ent_dim, tagbase->ncomps(), new_name);

      auto out =
          sync_array(ent_dim, as<Real>(tagbase)->array(), tagbase->ncomps());
      set_tag(ent_dim, new_name, out);

      change_tagToBoundary<Real> (ent_dim, tagbase->ncomps(), new_name);

      break;
    }
  }
}

void Mesh::change_all_bFieldsToMesh() {

  OMEGA_H_TIME_FUNCTION;
  for (Int ent_dim = 0; ent_dim < dim(); ++ent_dim) {
    for (Int t = 0; t < ntags(ent_dim); ++t) {
      auto tag = get_tag(ent_dim, t);
      auto const &name = tag->name();
      auto ncomps = tag->ncomps();

      if (is<I8>(tag)) {

        size_t found = name.find("_boundary");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToMesh<I8> (ent_dim, ncomps, name);
        }

      } else if (is<I32>(tag)) {

        size_t found = name.find("_boundary");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
          change_tagToMesh<I32> (ent_dim, ncomps, name);
        }

      } else if (is<I64>(tag)) {

        size_t found = name.find("_boundary");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToMesh<I64> (ent_dim, ncomps, name);
        }

      } else if (is<Real>(tag)) {

        size_t found = name.find("_boundary");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToMesh<Real> (ent_dim, ncomps, name);
        }
      }
    }
  }

  return;  
}

void Mesh::change_all_bFieldsToBoundary() {

  OMEGA_H_TIME_FUNCTION;
  for (Int ent_dim = 0; ent_dim < dim(); ++ent_dim) {
    for (Int t = 0; t < ntags(ent_dim); ++t) {
      auto tag = get_tag(ent_dim, t);
      auto const &name = tag->name();
      auto ncomps = tag->ncomps();

      if (is<I8>(tag)) {

        size_t found = name.find("_boundary");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToBoundary<I8> (ent_dim, ncomps, name);
        }

      } else if (is<I32>(tag)) {

        size_t found = name.find("_boundary");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
          change_tagToBoundary<I32> (ent_dim, ncomps, name);
        }

      } else if (is<I64>(tag)) {

        size_t found = name.find("_boundary");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToBoundary<I64> (ent_dim, ncomps, name);
        }

      } else if (is<Real>(tag)) {

        size_t found = name.find("_boundary");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToBoundary<Real> (ent_dim, ncomps, name);
        }
      }
    }
  }

  return;  
}

bool Mesh::has_allMeshTags() {

  OMEGA_H_TIME_FUNCTION;
  bool out = true;
  for (Int ent_dim = 0; ent_dim < dim(); ++ent_dim) {
    for (Int t = 0; t < ntags(ent_dim); ++t) {
      auto tag = get_tag(ent_dim, t);
      auto const &name = tag->name();
      size_t found = name.find("_boundary");
      if (found != std::string::npos) out = false;
    }
  }

  return out;

}

bool Mesh::has_anyBoundaryField() {

  return (!has_allMeshTags());

}

Adj Mesh::ask_revClass_downAdj (Int from, Int to) {
  auto rc = ask_revClass (from);
  auto ab2b = rc.ab2b;
  auto a2ab = rc.a2ab;
  auto n_gents = a2ab.size() - 1;
  auto nhighs = ab2b.size();
  auto down_ments = (ask_down(from, to)).ab2b;
  auto h2l_degree = element_degree (family(), from, to);

  Write<LO> g_hl2l (nhighs*h2l_degree);
  Write<LO> g2g_hl (n_gents + 1);

  auto f1 = OMEGA_H_LAMBDA(LO h) {
    LO h_id = ab2b[h];
    for (LO l = 0; l < h2l_degree; ++l) {
      g_hl2l[h*h2l_degree + l] = down_ments[h_id*h2l_degree + l];
    }
  };
  parallel_for(nhighs, f1, "createDownAb2b");

  auto f2 = OMEGA_H_LAMBDA(LO g) {
    g2g_hl[g] = a2ab[g]*h2l_degree;
  };
  parallel_for(n_gents + 1, f2, "createDownA2ab");

  return Adj(LOs(g2g_hl), LOs(g_hl2l));

}

}  // end namespace Omega_h
