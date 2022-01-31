#include <iostream>
#include "Omega_h_mesh.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_atomics.hpp"

namespace Omega_h {

#ifdef OMEGA_H_USE_CUDA
#endif

#define OMEGA_H_INST(T)                                                        \
  template void Mesh::change_tagTorc<T>(                                       \
      Int ent_dim, Int ncomps, std::string const& name, LOs class_ids);        \
  template void Mesh::change_tagToMesh<T>(                                     \
      Int ent_dim, Int ncomps, std::string const& name, LOs class_ids);        \
  template Read<T> Mesh::get_rcField_array<T>(                                 \
      Int dim, std::string const& name) const;                                 \
  template void Mesh::add_rcField<T>(                                          \
      Int dim, std::string const& name, Int ncomps);                           \
  template void Mesh::add_rcField<T>(                                          \
      LOs class_ids, Int dim, std::string const& name, Int ncomps);            \
  template void Mesh::add_rcField<T>(                                          \
      Int dim, std::string const& name, Int ncomps, Read<T> array);            \
  template void Mesh::set_rcField_array(                                       \
      Int dim, std::string const& name, Read<T> array); 
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

bool Mesh::has_revClass (Int edim) const {
  OMEGA_H_CHECK (has_ents(edim));
  return bool (revClass_[edim]);
}

Adj Mesh::get_revClass (Int edim) const {
  OMEGA_H_CHECK (has_ents(edim));
  OMEGA_H_CHECK (has_revClass (edim));
  return *(revClass_[edim]);
}

void sort_by_high_index(LOs const l2lh, Write<LO> const lh2h) {
  OMEGA_H_TIME_FUNCTION;
  LO const nl = l2lh.size() - 1;
  auto f = OMEGA_H_LAMBDA(LO const l) {
    LO const begin = l2lh[l];
    LO const end = l2lh[l + 1];
    for (LO j = begin; j < end; ++j) {
      LO k_min = j;
      GO min_h = lh2h[j];
      for (LO k = j + 1; k < end; ++k) {
        GO const h = lh2h[k];
        if (h < min_h) {
          k_min = k;
          min_h = h;
        }
      }
      swap2(lh2h[j], lh2h[k_min]);
    }
  };
  parallel_for(nl, std::move(f));
}

Adj Mesh::derive_revClass (Int edim, I8 should_sort) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK (has_ents(edim));

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

  if (should_sort > 0) {
    sort_by_high_index(a2ab_r, ab2b);
  }

  return Adj(a2ab_r, LOs(ab2b));

}

Adj Mesh::ask_revClass (Int edim) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK (has_ents(edim));
  if (has_revClass (edim)) {
    return get_revClass (edim);
  }
  Adj derived_rc = derive_revClass (edim);
  revClass_[edim] = std::make_shared<Adj>(derived_rc);
  return derived_rc;
}

Adj Mesh::ask_revClass (Int edim, LOs class_ids) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK (has_ents(edim));
  if (!class_ids.size()) {
    fprintf(stderr, "Model entity IDs cannot be empty\n");
    OMEGA_H_NORETURN(Adj());
  }
  auto edim_rc = ask_revClass(edim);
  auto n_gents = class_ids.size();
  auto max_gent_id = (edim_rc.a2ab.size() - 1);

  auto ab2b = edim_rc.ab2b;
  auto a2ab = edim_rc.a2ab;
  Write<LO> degree(n_gents, 0, "new_rc_degrees");

  auto count = OMEGA_H_LAMBDA (LO i) {
    auto gent = class_ids[i];
    if (gent <= max_gent_id) {
      auto start = a2ab[gent];
      auto end = a2ab[gent + 1];
      degree[i] = end - start;
    }
  };
  parallel_for(n_gents, std::move(count));

  auto total_ments = get_sum(LOs(degree));
  auto new_a2ab_r = offset_scan(LOs(degree), "new_rc_a2ab");
  Write<LO> new_ab2b(total_ments, 0, "new_rc_ab2b");

  auto get_values = OMEGA_H_LAMBDA (LO i) {
    auto gent = class_ids[i];
    if (gent <= max_gent_id) {
      auto start = a2ab[gent];
      auto end = a2ab[gent + 1];
      for (LO j = start; j < end; ++j) {
        new_ab2b[new_a2ab_r[i] + j - start] = ab2b[j];
      }
    }
  };
  parallel_for(n_gents, std::move(get_values));

  return Adj(new_a2ab_r, LOs(new_ab2b));
}

template <typename T>
void Mesh::change_tagTorc(Int ent_dim, Int ncomps, std::string const &name,
                          LOs class_ids) {

  OMEGA_H_TIME_FUNCTION;
  auto mesh_field = get_array<T>(ent_dim, name);
  OMEGA_H_CHECK (mesh_field.size() == nents(ent_dim)*ncomps);

  auto rc_ids = (ask_revClass(ent_dim)).ab2b;
  if (class_ids.exists()) rc_ids = (ask_revClass(ent_dim, class_ids)).ab2b;
  auto n_bEnts = rc_ids.size();
  Write<T> b_field(n_bEnts*ncomps);
  if ((ent_dim == 3) && (n_bEnts != nents(ent_dim))) {
    fprintf(stderr, "multiple model regions\n");
  }

  auto f = OMEGA_H_LAMBDA(LO i) {
    auto id = rc_ids[i];
    for (LO n = 0; n < ncomps; ++n) {
      b_field[i*ncomps + n] = mesh_field[id*ncomps + n];
    }
  };
  parallel_for(n_bEnts, f, "get_bdryField");

  set_tag<T>(ent_dim, name, Read<T>(b_field));

  return;
}

template <typename T>
void Mesh::change_tagToMesh(Int ent_dim, Int ncomps, std::string const &name,
                            LOs class_ids) {

  OMEGA_H_TIME_FUNCTION;
  auto rc_field = get_array<T>(ent_dim, name);
  auto n_ents = nents (ent_dim);
  OMEGA_H_CHECK (rc_field.size() <= n_ents*ncomps);

  auto rc_ids = (ask_revClass(ent_dim)).ab2b;
  if (class_ids.exists()) rc_ids = (ask_revClass(ent_dim, class_ids)).ab2b;
  auto n_bEnts = rc_ids.size();
  if ((ent_dim == 3) && (n_bEnts != n_ents)) {
    fprintf(stderr, "multiple model regions\n");
  }

  Write<T> mesh_field (n_ents*ncomps, OMEGA_H_INTERIOR_VAL);

  auto f = OMEGA_H_LAMBDA (LO i) {
    auto id = rc_ids[i];
    for (LO n = 0; n < ncomps; ++n) {
      if ((mesh_field[id*ncomps + n] - OMEGA_H_INTERIOR_VAL) < EPSILON) {
        mesh_field[id*ncomps + n] = rc_field[i*ncomps + n];
      }
    }
  };
  parallel_for(n_bEnts, f, "get_fieldFromBdry");

  set_tag<T>(ent_dim, name, Read<T>(mesh_field));

  return;
}

template <typename T>
Read<T> Mesh::get_rcField_array
  (Int ent_dim, std::string const& name) const {

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
    OMEGA_H_NORETURN(Read<T>());
  }

  std::string new_name = name;
  new_name.append("_rc");
  OMEGA_H_CHECK(has_tag(ent_dim, new_name));
  return get_tag<T>(ent_dim, new_name)->array();
}

template <typename T>
void Mesh::add_rcField(LOs class_ids, Int ent_dim, std::string const& name,
                       Int ncomps) {

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_rc");
  OMEGA_H_CHECK(!has_tag(ent_dim, new_name));
  check_dim2(ent_dim);
  check_tag_name(new_name);
  OMEGA_H_CHECK(ncomps >= 0);
  OMEGA_H_CHECK(ncomps <= Int(INT8_MAX));
  OMEGA_H_CHECK(tags_[ent_dim].size() < size_t(INT8_MAX));
  TagPtr ptr(new Tag<T>(new_name, ncomps, class_ids));
  tags_[ent_dim].push_back(std::move(ptr));

  return;
}

template <typename T>
void Mesh::add_rcField(Int ent_dim, std::string const& name, Int ncomps) {

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_rc");
  OMEGA_H_CHECK(!has_tag(ent_dim, new_name));
  add_tag<T>(ent_dim, new_name, ncomps);

  return;
}

template <typename T>
void Mesh::add_rcField(Int ent_dim, std::string const& name, Int ncomps,
                       Read<T> array) {

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_rc");
  OMEGA_H_CHECK(!has_tag(ent_dim, new_name));
  add_tag<T>(ent_dim, new_name, ncomps, array);

  return;
}

bool Mesh::has_rcField(Int ent_dim, std::string const& name) const {

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_rc");
  if (!has_ents(ent_dim)) return false;
  return has_tag(ent_dim, new_name);

}

void Mesh::remove_rcField(Int ent_dim, std::string const& name) {

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_rc");
  if (!has_ents(ent_dim)) return;
  remove_tag(ent_dim, new_name);

  return;
}

template <typename T>
void Mesh::set_rcField_array
  (Int ent_dim, std::string const& name, Read<T> array) {

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
  }

  std::string new_name = name;
  new_name.append("_rc");
  OMEGA_H_CHECK(has_tag(ent_dim, new_name));
  set_tag<T>(ent_dim, new_name, array);

  return;
}

void Mesh::reduce_rcField(Int ent_dim, std::string const& name,
     Omega_h_Op op) {

  std::string new_name = name;
  new_name.append("_rc");

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
  }
  else {
    OMEGA_H_CHECK(has_tag(ent_dim, new_name));
  }

  auto tagbase = get_tagbase(ent_dim, new_name);
  switch (tagbase->type()) {
    case OMEGA_H_I8: {

      if (nents(ent_dim)) {
        change_tagToMesh<I8> (ent_dim, tagbase->ncomps(), new_name,
                              tagbase->class_ids());
      }

      auto out = reduce_array(
          ent_dim, as<I8>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(ent_dim, new_name, out);

      change_tagTorc<I8> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());

      break;
    }
    case OMEGA_H_I32: {

      if (nents(ent_dim)) {
        change_tagToMesh<I32> (ent_dim, tagbase->ncomps(), new_name,
                              tagbase->class_ids());
      }

      auto out = reduce_array(
          ent_dim, as<I32>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(ent_dim, new_name, out);

      change_tagTorc<I32> (ent_dim, tagbase->ncomps(), new_name,
                           tagbase->class_ids());

      break;
    }
    case OMEGA_H_I64: {

      if (nents(ent_dim)) {
        change_tagToMesh<I64> (ent_dim, tagbase->ncomps(), new_name,
                               tagbase->class_ids());
      }

      auto out = reduce_array(
          ent_dim, as<I64>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(ent_dim, new_name, out);

      change_tagTorc<I64> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());

      break;
    }
    case OMEGA_H_F64: {

      if (nents(ent_dim)) { 
        change_tagToMesh<Real> (ent_dim, tagbase->ncomps(), new_name, 
                                tagbase->class_ids());
      }

      auto out = reduce_array(
          ent_dim, as<Real>(tagbase)->array(), tagbase->ncomps(), op);
      set_tag(ent_dim, new_name, out);

      change_tagTorc<Real> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());

      break;
    }
  }
}

void Mesh::sync_rcField(Int ent_dim, std::string const& name) {

  std::string new_name = name;
  new_name.append("_rc");

  size_t found = name.find("_rc");
  if (found != std::string::npos) {
    Omega_h_fail("duplicate suffix '_rc' at end of field name\n");
  }
  else {
    OMEGA_H_CHECK(has_tag(ent_dim, new_name));
  }

  auto tagbase = get_tagbase(ent_dim, new_name);
  switch (tagbase->type()) {
    case OMEGA_H_I8: {

      if (nents(ent_dim)) {
        change_tagToMesh<I8> (ent_dim, tagbase->ncomps(), new_name,
                              tagbase->class_ids());
      }

      auto out =
          sync_array(ent_dim, as<I8>(tagbase)->array(), tagbase->ncomps());
      set_tag(ent_dim, new_name, out);

      change_tagTorc<I8> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());

      break;
    }
    case OMEGA_H_I32: {

      if (nents(ent_dim)) { 
        change_tagToMesh<I32> (ent_dim, tagbase->ncomps(), new_name,
                               tagbase->class_ids());
      }

      auto out =
          sync_array(ent_dim, as<I32>(tagbase)->array(), tagbase->ncomps());
      set_tag(ent_dim, new_name, out);

      change_tagTorc<I32> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());

      break;
    }
    case OMEGA_H_I64: {

      if (nents(ent_dim)) { 
        change_tagToMesh<I64> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());
      }

      auto out =
          sync_array(ent_dim, as<I64>(tagbase)->array(), tagbase->ncomps());
      set_tag(ent_dim, new_name, out);

      change_tagTorc<I64> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());

      break;
    }
    case OMEGA_H_F64: {

      if (nents(ent_dim)) {
        change_tagToMesh<Real> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());
      }

      auto out =
          sync_array(ent_dim, as<Real>(tagbase)->array(), tagbase->ncomps());
      set_tag(ent_dim, new_name, out);

      change_tagTorc<Real> (ent_dim, tagbase->ncomps(), new_name,
                          tagbase->class_ids());

      break;
    }
  }
}

void Mesh::change_all_rcFieldsToMesh() {

  OMEGA_H_TIME_FUNCTION;
  for (Int ent_dim = 0; ent_dim < dim(); ++ent_dim) {
    for (Int t = 0; t < ntags(ent_dim); ++t) {
      auto tag = get_tag(ent_dim, t);
      auto const &name = tag->name();
      auto ncomps = tag->ncomps();
      auto class_ids = tag->class_ids();

      if (is<I8>(tag)) {

        size_t found = name.find("_rc");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToMesh<I8> (ent_dim, ncomps, name, class_ids);
        }

      } else if (is<I32>(tag)) {

        size_t found = name.find("_rc");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
          change_tagToMesh<I32> (ent_dim, ncomps, name, class_ids);
        }

      } else if (is<I64>(tag)) {

        size_t found = name.find("_rc");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToMesh<I64> (ent_dim, ncomps, name, class_ids);
        }

      } else if (is<Real>(tag)) {

        size_t found = name.find("_rc");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagToMesh<Real> (ent_dim, ncomps, name, class_ids);
        }
      }
    }
  }

  return;  
}

void Mesh::change_all_rcFieldsTorc() {

  OMEGA_H_TIME_FUNCTION;
  for (Int ent_dim = 0; ent_dim < dim(); ++ent_dim) {
    for (Int t = 0; t < ntags(ent_dim); ++t) {
      auto tag = get_tag(ent_dim, t);
      auto const &name = tag->name();
      auto ncomps = tag->ncomps();
      auto class_ids = tag->class_ids();

      if (is<I8>(tag)) {

        size_t found = name.find("_rc");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagTorc<I8> (ent_dim, ncomps, name, class_ids);
        }

      } else if (is<I32>(tag)) {

        size_t found = name.find("_rc");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
          change_tagTorc<I32> (ent_dim, ncomps, name, class_ids);
        }

      } else if (is<I64>(tag)) {

        size_t found = name.find("_rc");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagTorc<I64> (ent_dim, ncomps, name, class_ids);
        }

      } else if (is<Real>(tag)) {

        size_t found = name.find("_rc");
        if (found != std::string::npos) {
          if (nents(ent_dim)) 
            change_tagTorc<Real> (ent_dim, ncomps, name, class_ids);
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
      size_t found = name.find("_rc");
      if (found != std::string::npos) out = false;
    }
  }

  return out;

}

bool Mesh::has_anyrcField() {

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
