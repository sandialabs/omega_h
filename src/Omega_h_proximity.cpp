#include "Omega_h_confined.hpp"
#include "Omega_h_math.hpp"

template <Int dim>
static Reals get_edge_pad_isos(Mesh* mesh, Real max_size, Read<I8> edges_are_bridges) {
  auto coords = mesh->coords();
  auto edges2verts = mesh->ask_verts_of(EDGE);
  auto out = Write<Real>(mesh->nedges() * 2, max_size);
  auto f = LAMBDA(LO edge) {
    if (!edges_are_bridges[edge]) return;
    auto eev2v = gather_verts<2>(edges2verts, edge);
    auto eev2x = gather_vectors<2, dim>(coords, eev2v);
    auto l = norm(eev2x[1] - eev2x[0]);
    l = min2(l, max_size);
    out[edge * 2 + 0] = l;
    out[edge * 2 + 1] = l;
  };
  parallel_for(mesh->nedges(), f);
  return out;
}

template <Int dim>
static Reals get_tri_pad_isos(Mesh* mesh, Real max_size, Read<I8> edges_are_bridges) {
  auto coords = mesh->coords();
  auto tris2verts = mesh->ask_verts_of(TRI);
  auto tris2edges = mesh->ask_down(TRI, EDGE).ab2b;
  auto out = Write<Real>(mesh->ntris() * 3, max_size);
  auto f = LAMBDA(LO tri) {
    auto ttv2v = gather_verts<3>(tris2verts, tri);
    auto ttv2x = gather_vectors<3, dim>(coords, ttv2v);
    auto tte2e = gather_down<3>(tris2edges, tri);
    auto tte2b = gather_scalars<3>(edges_are_bridges, tte2e);
    for (Int ttv = 0; ttv < 3; ++ttv) {
      if (!(tte2b[(ttv + 0) % 3] && tte2b[(ttv + 2) % 3])) continue;
      auto o = ttv2x[ttv];
      auto a = ttv2x[(ttv + 1) % 3];
      auto b = ttv2x[(ttv + 2) % 3];
      auto oa = a - o;
      auto ab = b - a;
      auto nabsq = norm_squared(ab);
      auto proj = ab * ((ab * oa) / nabsq);
      auto d = oa - proj;
      auto lambda = ((ab * d) - (ab * oa)) / nabsq;
      if (!((0 <= lambda) && (lambda <= 1.0))) continue;
      auto h = norm(d);
      h = min2(h, max_size);
      out[tri * 3 + ttv] = h;
    }
  };
  parallel_for(mesh->nedges(), f);
  return out;
}

static Reals get_tet_pad_isos(Mesh* mesh, Real max_size, Read<I8> edges_are_bridges) {
  auto coords = mesh->coords();
  auto tets2verts = mesh->ask_verts_of(TET);
  auto tets2edges = mesh->ask_down(TET, EDGE).ab2b;
  auto out = Write<Real>(mesh->ntets() * 4, max_size);
  auto f = LAMBDA(LO tri) {
    auto ttv2v = gather_verts<4>(tris2verts, tri);
    auto ttv2x = gather_vectors<4, 3>(coords, ttv2v);
    auto tte2e = gather_down<6>(tris2edges, tri);
    auto tte2b = gather_scalars<6>(edges_are_bridges, tte2e);
    for (Int ttv = 0; ttv < 4; ++ttv) {
      Few<Int, 3> vve2tte;
      Few<Int, 3> vve2wd;
      for (Int vve = 0; vve < 3; ++vve) {
        vve2tte[vve] = UpTemplate<TET, VERT>::get(ttv, vve).up;
        vve2wd[vve] = UpTemplate<TET, VERT>::get(ttv, vve).which_down;
      }
      Few<I8, 3> vve2b;
      for (Int vve = 0; vve < 3; ++vve) vve2b[vve] = tte2b[vve2tte[vve]];
      if (minimum(vve2b) == 0) continue;
      auto o = ttv2x[ttv];
      Few<Int, 3> vve2ttv;
      for (Int vve = 0; vve < 3; ++vve) {
        vve2ttv[vve] = DownTemplate<TET, EDGE>(vve2tte[vve], 1 - vve2wd[vve]);
      }
      Few<Vector<3>, 3> vve2x;
      for (Int vve = 0; vve < 3; ++vve) vve2x[vve] = ttv2x[vve2ttv[vve]];
      auto n = normalize(cross(vve2x[1] - vve2x[0], vve2x[2] - vve2x[0]));
      auto oa = vve2x[0] - o;
      auto d = n * (n * oa);
    }
  };
  parallel_for(mesh->nedges(), f);
  return out;
}

static Reals get_pad_isos(Mesh* mesh, Real factor, LOs pads2elems) {
  if (mesh->dim() == 3) return get_pad_isos_dim<3>(mesh, factor, pads2elems);
  if (mesh->dim() == 2) return get_pad_isos_dim<2>(mesh, factor, pads2elems);
  NORETURN(Reals());
}

Reals get_proximity_isos(Mesh* mesh, Real factor, Real max_size) {
  CHECK(mesh->owners_have_all_upward(VERT));
  CHECK(mesh->owners_have_all_upward(EDGE));
  auto edges_are_bridges = find_bridge_edges(mesh);
  auto elems_are_pads = mark_up(mesh, EDGE, mesh->dim(), edges_are_bridges);
  auto elems_are_angle = find_angle_elems(mesh);
  auto elems_not_angle = invert_marks(elems_are_angle);
  elems_are_pads = land_each(elems_are_pads, elems_not_angle);
  auto pads2elems = collect_marked(elems_are_pads);
  auto pads2h = get_pad_isos(mesh, factor, pads2elems);
  auto elems2pads = invert_injective_map(pads2elems, mesh->nelems());
  auto verts2elems = mesh->ask_up(VERT, mesh->dim());
  auto vert_isos_w = Write<Real>(mesh->nverts());
  auto get_vert_values = LAMBDA(LO vert) {
    auto h = max_size;
    for (auto ve = verts2elems.a2ab[vert]; ve < verts2elems.a2ab[vert + 1];
         ++ve) {
      auto elem = verts2elems.ab2b[ve];
      auto pad = elems2pads[elem];
      if (pad < 0) continue;
      auto eh = pads2h[pad];
      h = min2(h, eh);
    }
    vert_isos_w[vert] = h;
  };
  parallel_for(mesh->nverts(), get_vert_values);
  auto vert_isos = Reals(vert_isos_w);
  return mesh->sync_array(VERT, vert_isos, 1);
}


