#include "surface.hpp"

#include "access.hpp"
#include "align.hpp"
#include "graph.hpp"
#include "loop.hpp"
#include "map.hpp"
#include "size.hpp"

namespace osh {

namespace surf {

namespace {

Reals get_triangle_normals(Mesh* mesh, LOs surf_tri2tri) {
  CHECK(mesh->dim() == 3);
  auto nsurf_tris = surf_tri2tri.size();
  auto fv2v = mesh->ask_verts_of(TRI);
  auto coords = mesh->coords();
  Write<Real> normals(nsurf_tris * 3);
  auto lambda = LAMBDA(LO surf_tri) {
    auto f = surf_tri2tri[surf_tri];
    auto v = gather_verts<3>(fv2v, f);
    auto x = gather_vectors<3, 3>(coords, v);
    auto b = simplex_basis<3, 2>(x);
    auto n = normalize(cross(b[0], b[1]));
    set_vector(normals, surf_tri, n);
  };
  parallel_for(nsurf_tris, lambda);
  return normals;
}

Reals get_edge_normals(Mesh* mesh, LOs surf_edge2edge) {
  CHECK(mesh->dim() == 2);
  auto nsurf_edges = surf_edge2edge.size();
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto coords = mesh->coords();
  Write<Real> normals(nsurf_edges * 2);
  auto lambda = LAMBDA(LO surf_edge) {
    auto e = surf_edge2edge[surf_edge];
    auto v = gather_verts<2>(ev2v, e);
    auto x = gather_vectors<2, 2>(coords, v);
    auto b = simplex_basis<2, 1>(x);
    auto n = normalize(perp(b[0]));
    set_vector(normals, surf_edge, n);
  };
  parallel_for(nsurf_edges, lambda);
  return normals;
}

template <Int dim>
Reals get_hinge_angles_tmpl(Mesh* mesh, Reals surf_side_normals,
                            LOs surf_hinge2hinge, LOs side2surf_side) {
  auto nsurf_hinges = surf_hinge2hinge.size();
  auto hinges2sides = mesh->ask_up(dim - 2, dim - 1);
  auto hinges2hinge_sides = hinges2sides.a2ab;
  auto hinge_sides2sides = hinges2sides.ab2b;
  Write<Real> angles(nsurf_hinges);
  auto f = LAMBDA(LO surf_hinge) {
    auto hinge = surf_hinge2hinge[surf_hinge];
    auto begin = hinges2hinge_sides[hinge];
    auto end = hinges2hinge_sides[hinge + 1];
    Int i = 0;
    Vector<dim> n[2];
    for (auto hs = begin; hs < end; ++hs) {
      auto s = hinge_sides2sides[hs];
      auto ss = side2surf_side[s];
      if (-1 == ss) continue;
      n[i++] = get_vector<dim>(surf_side_normals, ss);
    }
    angles[surf_hinge] = acos(n[0] * n[1]);
  };
  parallel_for(nsurf_hinges, f);
  return angles;
}

}  // end anonymous namespace

Reals get_side_normals(Mesh* mesh, LOs surf_side2side) {
  if (mesh->dim() == 3) {
    return get_triangle_normals(mesh, surf_side2side);
  } else {
    return get_edge_normals(mesh, surf_side2side);
  }
}

Reals get_hinge_angles(Mesh* mesh, Reals surf_side_normals,
                       LOs surf_hinge2hinge, LOs side2surf_side) {
  if (mesh->dim() == 3) {
    return get_hinge_angles_tmpl<3>(mesh, surf_side_normals, surf_hinge2hinge,
                                    side2surf_side);
  } else {
    return get_hinge_angles_tmpl<2>(mesh, surf_side_normals, surf_hinge2hinge,
                                    side2surf_side);
  }
}

/*
 * Max, Nelson.
 * "Weights for computing vertex normals from facet normals."
 * Journal of Graphics Tools 4.2 (1999): 1-6.
 */

static Reals tri_vert_normal_weights(Mesh* mesh, LOs surf_tri2tri) {
  auto ntris = mesh->nents(TRI);
  auto tri2surf_tri = invert_injective_map(surf_tri2tri, ntris);
  auto v2f = mesh->ask_up(VERT, TRI);
  auto fv2v = mesh->ask_verts_of(TRI);
  auto coords = mesh->coords();
  auto v2vf = v2f.a2ab;
  auto nvf = v2vf.last();
  auto vf2f = v2f.ab2b;
  auto weights = Write<Real>(nvf);
  auto func = LAMBDA(LO v) {
    for (auto vf = v2vf[v]; vf < v2vf[v + 1]; ++vf) {
      auto f = vf2f[vf];
      if (tri2surf_tri[f] < 0) {
        weights[vf] = 0.0;
        continue;
      }
      auto code = v2f.codes[vf];
      auto ffv = code_which_down(code);
      auto rot = rotation_to_first<3>(ffv);
      auto ffv2v_0 = gather_verts<3>(fv2v, f);
      Few<LO, 3> ffv2v;
      rotate_adj<3>(rot, &ffv2v_0[0], &ffv2v[0]);
      auto ffv2p = gather_vectors<3, 3>(coords, ffv2v);
      auto b = simplex_basis<3, 2>(ffv2p);
      auto w =
          norm(cross(b[0], b[1])) / (norm_squared(b[0]) * norm_squared(b[1]));
      weights[vf] = w;
    }
  };
  parallel_for(mesh->nverts(), func);
  return weights;
}

/* following the same derivation as Max did above,
 * we find that the weights for edges should be
 * one over the length of the edge.
 */

static Reals edge_vert_normal_weights(Mesh* mesh, LOs surf_edge2edge) {
  auto surf_edge2len = measure_edges_real(mesh, surf_edge2edge);
  auto nsurf_edges = surf_edge2edge.size();
  auto surf_edge2weight_w = Write<Real>(nsurf_edges);
  auto f = LAMBDA(LO se) { surf_edge2weight_w[se] = 1.0 / surf_edge2len[se]; };
  parallel_for(nsurf_edges, f);
  auto surf_edge2weight = Reals(surf_edge2weight_w);
  auto edge2weight =
      map_onto(surf_edge2weight, surf_edge2edge, mesh->nedges(), 0.0, 1);
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto ve2e = v2e.ab2b;
  auto ve2weight = unmap(ve2e, edge2weight, 1);
  return ve2weight;
}

static Reals side_vert_normal_weights(Mesh* mesh, LOs surf_side2side) {
  if (mesh->dim() == 3) return tri_vert_normal_weights(mesh, surf_side2side);
  if (mesh->dim() == 2) return edge_vert_normal_weights(mesh, surf_side2side);
  NORETURN(Reals());
}

Reals get_vert_normals(Mesh* mesh, LOs surf_side2side, Reals surf_side_normals,
                       LOs surf_vert2vert) {
  CHECK(mesh->owners_have_all_upward(VERT));
  auto weights = side_vert_normal_weights(mesh, surf_side2side);
  auto dim = mesh->dim();
  auto v2s = mesh->ask_graph(VERT, dim - 1);
  auto nsides = mesh->nents(dim - 1);
  auto side_normals =
      map_onto(surf_side_normals, surf_side2side, nsides, 0.0, dim);
  auto vert_normals = graph_weighted_average(v2s, weights, side_normals, dim);
  vert_normals = mesh->sync_array(VERT, vert_normals, dim);
  auto surf_vert_normals = unmap(surf_vert2vert, vert_normals, dim);
  return normalize_vectors(surf_vert_normals, dim);
}

}  // end namespace surf

}  // end namespace osh
