#include "Omega_h_transfer.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_motion.hpp"
#include "Omega_h_shape.hpp"

namespace Omega_h {

static bool should_transfer_motion_linear(
    Mesh* mesh, TransferOpts const& opts, TagBase const* tb) {
  if (tb->name() == "warp") return false;
  return should_interpolate(mesh, opts, VERT, tb) ||
         is_metric(mesh, opts, VERT, tb);
}

LinearPack pack_linearized_fields(Mesh* mesh, TransferOpts const& opts) {
  Int ncomps = 0;
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    auto tb = mesh->get_tag(VERT, i);
    if (should_transfer_motion_linear(mesh, opts, tb)) {
      ncomps += tb->ncomps();
    }
  }
  auto out_w = Write<Real>(mesh->nverts() * ncomps);
  Int offset = 0;
  Int metric_offset = -1;
  Int coords_offset = -1;
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    auto tb = mesh->get_tag(VERT, i);
    if (!should_transfer_motion_linear(mesh, opts, tb)) continue;
    auto t = dynamic_cast<Tag<Real> const*>(tb);
    auto in = t->array();
    if (is_metric(mesh, opts, VERT, tb)) {
      in = linearize_metrics(mesh->nverts(), in);
    }
    auto ncomps_in = tb->ncomps();
    auto f = OMEGA_H_LAMBDA(LO v) {
      for (Int c = 0; c < ncomps_in; ++c) {
        out_w[v * ncomps + offset + c] = in[v * ncomps_in + c];
      }
    };
    parallel_for(mesh->nverts(), f);
    if (tb->name() == "metric") metric_offset = offset;
    if (tb->name() == "coordinates") coords_offset = offset;
    offset += ncomps_in;
  }
  return {out_w, ncomps, metric_offset, coords_offset};
}

void unpack_linearized_fields(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, Reals data, Read<I8> verts_are_keys) {
  OMEGA_H_CHECK(data.size() % new_mesh->nverts() == 0);
  auto ncomps = data.size() / new_mesh->nverts();
  Int offset = 0;
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tb = old_mesh->get_tag(VERT, i);
    if (!should_transfer_motion_linear(old_mesh, opts, tb)) continue;
    auto ncomps_out = tb->ncomps();
    auto out_w = Write<Real>(new_mesh->nverts() * ncomps_out);
    auto f = OMEGA_H_LAMBDA(LO v) {
      for (Int c = 0; c < ncomps_out; ++c) {
        out_w[v * ncomps_out + c] = data[v * ncomps + offset + c];
      }
    };
    parallel_for(new_mesh->nverts(), f);
    auto out = Reals(out_w);
    if (is_metric(old_mesh, opts, VERT, tb)) {
      out = delinearize_metrics(old_mesh->nverts(), out);
    }
    auto t = dynamic_cast<Tag<Real> const*>(tb);
    auto prev = t->array();
    out_w = deep_copy(prev);
    auto f2 = OMEGA_H_LAMBDA(LO v) {
      if (!verts_are_keys[v]) return;
      for (Int c = 0; c < ncomps_out; ++c) {
        out_w[v * ncomps_out + c] = out[v * ncomps_out + c];
      }
    };
    parallel_for(new_mesh->nverts(), f2);
    out = Reals(out_w);
    if (new_mesh->has_tag(VERT, tb->name())) {
      new_mesh->set_tag(VERT, tb->name(), out);
    } else {
      new_mesh->add_tag(VERT, tb->name(), ncomps_out, out);
    }
    offset += ncomps_out;
  }
  OMEGA_H_CHECK(offset == ncomps);
}

}  // end namespace Omega_h
