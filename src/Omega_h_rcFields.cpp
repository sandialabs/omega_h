#include <algorithm>
#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_atomics.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_timer.hpp"

namespace Omega_h {

#ifdef OMEGA_H_USE_CUDA
#endif

#define OMEGA_H_INST(T)                                                        \
  template void Mesh::change_tagTorc<T>(                                       \
      Int ent_dim, Int ncomps, std::string const& name, LOs class_ids, bool);  \
  template void Mesh::change_tagToMesh<T>(                                     \
      Int ent_dim, Int ncomps, std::string const& name, LOs class_ids, bool);  \
  template Read<T> Mesh::get_rcField_array<T>(                                 \
      Int dim, std::string const& name) const;                                 \
  template void Mesh::add_rcField<T>(                                          \
      Int dim, std::string const& name, Int ncomps);                           \
  template void Mesh::add_rcField<T>(                                          \
      LOs class_ids, Int dim, std::string const& name, Int ncomps);            \
  template void Mesh::add_rcField<T>(                                          \
      Int dim, std::string const& name, Int ncomps, Read<T> array);            \
  template void Mesh::set_rcField_array(                                       \
      Int dim, std::string const& name, Read<T> array);                        \
  template Read<T> Mesh::get_rc_mesh_array(                                    \
      Int ent_dim, Int ncomps, std::string const& name, LOs class_ids);        \
  template Read<T> Mesh::get_rc_array_from_mesh_array(Int ent_dim, Int ncomps, \
      std::string const& name, LOs class_ids, Read<T> mesh_array);             \
  template void Mesh::set_rc_from_mesh_array(Int ent_dim, Int ncomps,          \
      LOs class_ids, std::string const& name, Read<T> array);                  \
  template std::unique_ptr<Tag<T>> Mesh::get_rc_mesh_tag_from_rc_tag(          \
      Int, Tag<T> const*);                                                     \
  template Read<T> Mesh::get_rc_array(Int dim, std::string const& name) const; \
  template Read<T> Mesh::get_rc_mesh_array_from_rc_array(                      \
      Int ent_dim, Int ncomps, LOs class_ids, Read<T> rc_field);
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

static std::string get_rc_name(std::string name) {
  size_t found = name.find("_rc");
  if (found == std::string::npos) {
    return name.append("_rc");
  }
  return name;
}
template <typename T>
Read<T> Mesh::get_rc_array(Int dim, std::string const& name) const {
  auto tag_itr = rc_tag_iter(dim, name);
  OMEGA_H_CHECK(tag_itr.had_tag);
  return as<T>(tag_itr.it->get())->array();
}

Int Mesh::nrctags(Int dim) const { return rc_field_tags_[dim].size(); }

bool Mesh::has_revClass(Int edim) const {
  OMEGA_H_CHECK(has_ents(edim));
  return bool(revClass_[edim]);
}

Adj Mesh::get_revClass(Int edim) const {
  OMEGA_H_CHECK(has_ents(edim));
  OMEGA_H_CHECK(has_revClass(edim));
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

Adj Mesh::derive_revClass(Int edim, I8 should_sort) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(has_ents(edim));

  auto class_ids_all = get_array<ClassId>(edim, "class_id");
  auto class_dim = get_array<I8>(edim, "class_dim");
  // copy_if
  Write<LO> class_ids_w(nents(edim), -1, "edim_classIds");
  auto edim_classid = OMEGA_H_LAMBDA(LO i) {
    if (class_dim[i] == edim) {
      class_ids_w[i] = class_ids_all[i];
    }
  };
  parallel_for(nents(edim), std::move(edim_classid));
  auto class_ids = LOs(class_ids_w);

  auto const n_gents = get_max(class_ids) + 1;
  Write<LO> degree(n_gents, 0, "rc_degrees");
  auto count_degree = OMEGA_H_LAMBDA(LO i) {
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

  auto get_values = OMEGA_H_LAMBDA(LO i) {
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

  return {a2ab_r, LOs(ab2b)};
}

Adj Mesh::ask_revClass(Int edim) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(has_ents(edim));
  if (has_revClass(edim)) {
    return get_revClass(edim);
  }
  Adj derived_rc = derive_revClass(edim);
  revClass_[edim] = std::make_shared<Adj>(derived_rc);
  return derived_rc;
}

Adj Mesh::ask_revClass(Int edim, LOs class_ids) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(has_ents(edim));
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

  auto count = OMEGA_H_LAMBDA(LO i) {
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

  auto get_values = OMEGA_H_LAMBDA(LO i) {
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
  return {new_a2ab_r, LOs(new_ab2b)};
}

template <typename T>
Read<T> Mesh::get_rc_array_from_mesh_array(Int ent_dim, Int ncomps,
    std::string const& name, LOs class_ids, Read<T> mesh_field) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(mesh_field.size() == nents(ent_dim) * ncomps);
  auto rc_ids = (ask_revClass(ent_dim)).ab2b;
  if (class_ids.exists()) rc_ids = (ask_revClass(ent_dim, class_ids)).ab2b;
  auto n_bEnts = rc_ids.size();
  Write<T> b_field(n_bEnts * ncomps);
  if ((ent_dim == 3) && (n_bEnts != nents(ent_dim))) {
    fprintf(stderr, "multiple model regions\n");
  }

  auto f = OMEGA_H_LAMBDA(LO i) {
    auto id = rc_ids[i];
    for (LO n = 0; n < ncomps; ++n) {
      b_field[i * ncomps + n] = mesh_field[id * ncomps + n];
    }
  };
  parallel_for(n_bEnts, f, "get_bdryField");
  return b_field;
}

void Mesh::set_rc_from_mesh_array(Int ent_dim, Int ncomps, LOs class_ids,
    std::string const& name, Read<T> array) {
  OMEGA_H_TIME_FUNCTION;
  auto tag_itr = rc_tag_iter(ent_dim,name);
  auto b_field =
      get_rc_array_from_mesh_array(ent_dim, ncomps, name, class_ids, array);
  if(!tag_itr.had_tag) {
    add_rcField<T>(class_ids,ent_dim,name,ncomps);
  }
  set_rcField_array(ent_dim, name, b_field);
}

// TODO take tag as parameter, return pointer to new tag
template <typename T>
void Mesh::change_tagTorc(Int ent_dim, Int ncomps, std::string const& name,
    LOs class_ids, bool remove) {
  OMEGA_H_TIME_FUNCTION;
  auto rc_name = get_rc_name(name);
  auto mesh_field = get_array<T>(ent_dim, rc_name);
  OMEGA_H_CHECK(mesh_field.size() == nents(ent_dim) * ncomps);
  auto b_field = get_rc_array_from_mesh_array(
      ent_dim, ncomps, rc_name, class_ids, mesh_field);
  add_rcField<T>(class_ids, ent_dim, rc_name, ncomps);
  set_rcField_array(ent_dim, rc_name, b_field);
  if (remove) {
    remove_tag(ent_dim, rc_name);
  }
}

// TODO make const after ask_revClass is made const
template <typename T>
Read<T> Mesh::get_rc_mesh_array_from_rc_array(
    Int ent_dim, Int ncomps, LOs class_ids, Read<T> rc_field) {
  OMEGA_H_TIME_FUNCTION;
  auto n_ents = nents(ent_dim);
  // if there are no entities in the requested dimension return an empty array
  if (n_ents == 0) {
    return {};
  }

  OMEGA_H_CHECK(rc_field.size() <= n_ents * ncomps);

  auto rc_ids = class_ids.exists() ? (ask_revClass(ent_dim, class_ids)).ab2b
                                   : (ask_revClass(ent_dim)).ab2b;
  auto class_id_size = class_ids.exists() ? class_ids.size() : -1;
  auto n_bEnts = rc_ids.size();
  OMEGA_H_CHECK(rc_field.size() == n_bEnts * ncomps);
  if ((ent_dim == 3) && (n_bEnts != n_ents)) {
    fprintf(stderr, "multiple model regions\n");
  }

  Write<T> mesh_field(n_ents * ncomps, OMEGA_H_INTERIOR_VAL);

  auto f = OMEGA_H_LAMBDA(LO i) {
    auto id = rc_ids[i];
    for (LO n = 0; n < ncomps; ++n) {
      if ((mesh_field[id * ncomps + n] - OMEGA_H_INTERIOR_VAL) < EPSILON) {
        mesh_field[id * ncomps + n] = rc_field[i * ncomps + n];
      }
    }
  };
  parallel_for(n_bEnts, f, "get_fieldFromBdry");
  return Read<T>(mesh_field);
}

template <typename T>
Read<T> Mesh::get_rc_mesh_array(
    Int ent_dim, Int ncomps, std::string const& name, LOs class_ids) {
  OMEGA_H_TIME_FUNCTION;
  auto rc_field = get_rc_array<T>(ent_dim, name);
  return get_rc_mesh_array_from_rc_array(ent_dim, ncomps, class_ids, rc_field);
}

// TODO take tag as parameter, return pointer to new tag
template <typename T>
void Mesh::change_tagToMesh(Int ent_dim, Int ncomps, std::string const& name,
    LOs class_ids, bool remove) {
  auto mesh_field = get_rc_mesh_array<T>(ent_dim, ncomps, name, class_ids);
  add_tag<T>(ent_dim, name, ncomps, Read<T>(mesh_field));
  if (remove) {
    remove_rcField(ent_dim, name);
  }
}

template <typename T>
Read<T> Mesh::get_rcField_array(Int ent_dim, std::string const& name) const {
  auto new_name = get_rc_name(name);
  auto tag_itr = rc_tag_iter(ent_dim, new_name);
  OMEGA_H_CHECK(tag_itr.had_tag == true);
  return as<T>(tag_itr.it->get())->array();
}

// TODO make const once get_rc_mesh_array_from_rc_arraay made const
template <typename T>
std::unique_ptr<Tag<T>> Mesh::get_rc_mesh_tag_from_rc_tag(
    Int ent_dim, Tag<T> const* tag) {
  const auto& name = tag->name();
  const auto ncomps = tag->ncomps();
  const auto class_ids = tag->class_ids();
  auto new_name = get_rc_name(name);
  OMEGA_H_CHECK(ncomps >= 0);
  OMEGA_H_CHECK(ncomps <= Int(INT8_MAX));
  auto new_tag = std::make_unique<Tag<T>>(new_name, ncomps, class_ids);
  auto rc_array = tag->array();
  auto mesh_array =
      get_rc_mesh_array_from_rc_array(ent_dim, ncomps, class_ids, rc_array);
  new_tag->set_array(mesh_array);
  return new_tag;
}

std::unique_ptr<TagBase> Mesh::get_rc_mesh_tag_from_rc_tag(
    Int ent_dim, TagBase const* tag) {
  auto new_tag = apply_to_omega_h_types(tag->type(), [&](auto t) {
    using T = decltype(t);
    return std::unique_ptr<TagBase>{
        get_rc_mesh_tag_from_rc_tag(ent_dim, as<T>(tag))};
  });
  return new_tag;
}

// TODO remove name argument and use name from TagpPtr
void Mesh::add_rcField(Int ent_dim, std::string const& name, TagPtr tag) {
  auto new_name = get_rc_name(name);
  OMEGA_H_CHECK(tag->name() == new_name);
  auto tag_itr = rc_tag_iter(ent_dim, new_name);
  check_dim2(ent_dim);
  check_tag_name(new_name);
  OMEGA_H_CHECK(rc_field_tags_[ent_dim].size() < size_t(INT8_MAX));
  if (tag_itr.had_tag) {
    *tag_itr.it = std::move(tag);
  } else {
    rc_field_tags_[ent_dim].push_back(std::move(tag));
  }
}

template <typename T>
void Mesh::add_rcField(
    LOs class_ids, Int ent_dim, std::string const& name, Int ncomps) {
  auto new_name = get_rc_name(name);
  OMEGA_H_CHECK(ncomps >= 0);
  OMEGA_H_CHECK(ncomps <= Int(INT8_MAX));
  auto tag = std::make_shared<Tag<T>>(new_name, ncomps, class_ids);
  add_rcField(ent_dim, new_name, std::move(tag));
}

template <typename T>
void Mesh::add_rcField(Int ent_dim, std::string const& name, Int ncomps) {
  auto new_name = get_rc_name(name);
  auto tag = std::make_shared<Tag<T>>(new_name, ncomps);
  add_rcField(ent_dim, new_name, std::move(tag));
}

template <typename T>
void Mesh::add_rcField(
    Int ent_dim, std::string const& name, Int ncomps, Read<T> array) {
  auto new_name = get_rc_name(name);
  // FIXME should have version without class_ids array...
  auto tag = std::make_shared<Tag<T>>(new_name, ncomps);
  tag->set_array(array);
  add_rcField(ent_dim, new_name, std::move(tag));
}

bool Mesh::has_rcField(Int ent_dim, std::string const& name) const {
  auto new_name = get_rc_name(name);
  auto tag_itr = rc_tag_iter(ent_dim, new_name);
  return tag_itr.had_tag;
}

void Mesh::remove_rcField(Int ent_dim, std::string const& name) {
  auto new_name = get_rc_name(name);
  auto tag_itr= rc_tag_iter(ent_dim, new_name);
  if (tag_itr.had_tag) {
    rc_field_tags_[ent_dim].erase(tag_itr.it);
  }
}

template <typename T>
void Mesh::set_rcField_array(
    Int ent_dim, std::string const& name, Read<T> array) {
  auto new_name = get_rc_name(name);
  auto tag_itr = rc_tag_iter(ent_dim, new_name);
  OMEGA_H_CHECK(tag_itr.had_tag);
  auto class_ids = (*tag_itr.it)->class_ids();
  // all rc_tags should have
  const auto ncomps = (*tag_itr.it)->ncomps();
  auto tag = std::make_shared<Tag<T>>(new_name, ncomps, class_ids);
  tag->set_array(array);
  add_rcField(ent_dim, new_name, std::move(tag));
}

void Mesh::reduce_rcField(Int ent_dim, std::string const& name, Omega_h_Op op) {
  auto new_name = get_rc_name(name);

  auto tag_itr= rc_tag_iter(ent_dim, new_name);
  OMEGA_H_CHECK(tag_itr.had_tag);
  // auto tagbase
  auto* tagbase = tag_itr.it->get();
  const auto class_ids = tagbase->class_ids();
  const auto ncomps = tagbase->ncomps();

  auto f = [&](auto t) {
    using T = decltype(t);
    auto mesh_array =
        get_rc_mesh_array<T>(ent_dim, ncomps, new_name, class_ids);
    auto out = reduce_array(ent_dim, mesh_array, ncomps, op);
    set_rc_from_mesh_array(ent_dim, ncomps, class_ids, new_name, out);
  };
  apply_to_omega_h_types(tagbase->type(), f);
}

void Mesh::sync_rcField(Int ent_dim, std::string const& name) {
  auto new_name = get_rc_name(name);
  auto tag_itr = rc_tag_iter(ent_dim, new_name);
  OMEGA_H_CHECK(tag_itr.had_tag);

  const auto ncomps = (*tag_itr.it)->ncomps();
  const auto class_ids = (*tag_itr.it)->class_ids();
  auto f = [&](auto t) {
    using T = decltype(t);
    auto mesh_array =
        get_rc_mesh_array<T>(ent_dim, ncomps, new_name, class_ids);
    auto out = sync_array(ent_dim, mesh_array, ncomps);
    set_rc_from_mesh_array(ent_dim, ncomps, class_ids, new_name, out);
  };
  apply_to_omega_h_types((*tag_itr.it)->type(), f);
}

bool Mesh::change_all_rcFieldsToMesh() {
  OMEGA_H_TIME_FUNCTION;
  bool changed = false;
  for (Int ent_dim = 0; ent_dim <= dim(); ++ent_dim) {
    for (const auto& rc_tag : rc_field_tags_[ent_dim]) {
      changed = true;
      OMEGA_H_CHECK(rc_tag != nullptr);
      auto const& name = rc_tag->name();
      auto ncomps = rc_tag->ncomps();
      auto class_ids = rc_tag->class_ids();
      auto f = [&](auto t) {
        using T = decltype(t);
        // TODO can change_tagToMesh just take a tag? because we get all the
        // data from the tag anyways...
        change_tagToMesh<T>(ent_dim, ncomps, name, class_ids, false);
      };
      apply_to_omega_h_types(rc_tag->type(), f);
      OMEGA_H_CHECK(has_tag(ent_dim, name));
    }
    rc_field_tags_[ent_dim].clear();
    OMEGA_H_CHECK(rc_field_tags_[ent_dim].size() == 0);
  }
  return changed;
}

bool Mesh::change_all_rcFieldsTorc() {
  OMEGA_H_TIME_FUNCTION;
  bool changed = false;
  for (Int ent_dim = 0; ent_dim <= dim(); ++ent_dim) {
    for (const auto& tag : tags_[ent_dim]) {
      if (is_rc_tag(tag->name())) {
        changed = true;
        OMEGA_H_CHECK(tag != nullptr);
        auto const& name = tag->name();
        auto ncomps = tag->ncomps();
        auto class_ids = tag->class_ids();
        auto f = [&](auto t) {
          using T = decltype(t);
          // TODO can change_tagToMesh just take a tag? because we get all the
          // data from the tag anyways...
          change_tagTorc<T>(ent_dim, ncomps, name, class_ids, false);
        };
        apply_to_omega_h_types(tag->type(), f);
      }
    }
    tags_[ent_dim].erase(
        std::remove_if(tags_[ent_dim].begin(), tags_[ent_dim].end(),
            [](const auto& tag) {
              return is_rc_tag(tag->name());
            }),
        tags_[ent_dim].end());
  }
  return changed;
}

Adj Mesh::ask_revClass_downAdj(Int from, Int to) {
  auto rc = ask_revClass(from);
  auto ab2b = rc.ab2b;
  auto a2ab = rc.a2ab;
  auto n_gents = a2ab.size() - 1;
  auto nhighs = ab2b.size();
  auto down_ments = (ask_down(from, to)).ab2b;
  auto h2l_degree = element_degree(family(), from, to);

  Write<LO> g_hl2l(nhighs * h2l_degree);
  Write<LO> g2g_hl(n_gents + 1);

  auto f1 = OMEGA_H_LAMBDA(LO h) {
    LO h_id = ab2b[h];
    for (LO l = 0; l < h2l_degree; ++l) {
      g_hl2l[h * h2l_degree + l] = down_ments[h_id * h2l_degree + l];
    }
  };
  parallel_for(nhighs, f1, "createDownAb2b");

  auto f2 = OMEGA_H_LAMBDA(LO g) { g2g_hl[g] = a2ab[g] * h2l_degree; };
  parallel_for(n_gents + 1, f2, "createDownA2ab");
  return {LOs(g2g_hl), LOs(g_hl2l)};
}
bool is_rc_tag(std::string const& name) {
  return (name.find("_rc") != std::string::npos);
}

}  // end namespace Omega_h
