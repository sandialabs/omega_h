#include "Omega_h_motion.hpp"
#include "access.hpp"
#include "loop.hpp"
#include "array.hpp"
#include "metric.hpp"
#include "size.hpp"

namespace Omega_h {

static bool should_transfer(TagBase const* tb) {
  return (tb->xfer() == OMEGA_H_LINEAR_INTERP && tb->name() != "warp") ||
    tb->xfer() == OMEGA_H_METRIC || tb->xfer() == OMEGA_H_SIZE;
}

LinearPack pack_linearized_fields(Mesh* mesh) {
  Int ncomps = 0;
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    auto tb = mesh->get_tag(VERT, i);
    if (should_transfer(tb)) {
      ncomps += tb->ncomps();
    }
  }
  auto out_w = Write<Real>(mesh->nverts() * ncomps);
  Int offset = 0;
  Int metric_offset = -1;
  Int coords_offset = -1;
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    auto tb = mesh->get_tag(VERT, i);
    if (!should_transfer(tb)) continue;
    auto t = dynamic_cast<Tag<Real> const*>(tb);
    auto in = t->array();
    if (tb->xfer() == OMEGA_H_METRIC) in = linearize_metrics(mesh->dim(), in);
    else if (tb->xfer() == OMEGA_H_SIZE) in = linearize_isos(in);
    auto ncomps_in = tb->ncomps();
    auto f = LAMBDA(LO v) {
      for (Int i = 0; i < ncomps_in; ++i) {
        out_w[v * ncomps + offset + i] = in[v * ncomps_in + i];
      }
    };
    parallel_for(mesh->nverts(), f);
    if (tb->name() == "metric" || tb->name() == "size") metric_offset = offset;
    if (tb->name() == "coordinates") coords_offset = offset;
    offset += ncomps_in;
  }
  return { out_w, ncomps, metric_offset, coords_offset };
}

void unpack_linearized_fields(Mesh* old_mesh, Mesh* new_mesh, Reals data,
    Read<I8> verts_are_keys) {
  CHECK(data.size() % new_mesh->nverts() == 0);
  auto ncomps = data.size() / new_mesh->nverts();
  Int offset = 0;
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tb = old_mesh->get_tag(VERT, i);
    if (!should_transfer(tb)) continue;
    auto ncomps_out = tb->ncomps();
    auto out_w = Write<Real>(new_mesh->nverts() * ncomps_out);
    auto f = LAMBDA(LO v) {
      for (Int i = 0; i < ncomps_out; ++i) {
        out_w[v * ncomps_out + i] = data[v * ncomps + offset + i];
      }
    };
    parallel_for(new_mesh->nverts(), f);
    auto out = Reals(out_w);
    if (tb->xfer() == OMEGA_H_METRIC) out = delinearize_metrics(old_mesh->dim(), out);
    else if (tb->xfer() == OMEGA_H_SIZE) out = delinearize_isos(out);
    auto t = dynamic_cast<Tag<Real> const*>(tb);
    auto prev = t->array();
    out_w = deep_copy(prev);
    auto f2 = LAMBDA(LO v) {
      if (!verts_are_keys[v]) return;
      for (Int i = 0; i < ncomps_out; ++i) {
        out_w[v * ncomps_out + i] = out[v * ncomps_out + i];
      }
    };
    parallel_for(new_mesh->nverts(), f2);
    out = Reals(out_w);
    if (new_mesh->has_tag(VERT, tb->name())) {
      new_mesh->set_tag(VERT, tb->name(), out);
    } else {
      new_mesh->add_tag(
          VERT, tb->name(), ncomps_out, tb->xfer(), tb->outflags(), out);
    }
    offset += ncomps_out;
  }
  CHECK(offset == ncomps);
}

} // end namespace Omega_h
