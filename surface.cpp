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
    /* right now omega_h generates edges counter-clockwise
     * around triangles, so the outward-facing normal is the
     * negative of what perp() gives
     */
    auto n = -normalize(perp(b[0]));
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
  switch (mesh->dim()) {
    case 3:
      return get_triangle_normals(mesh, surf_side2side);
    case 2:
      return get_edge_normals(mesh, surf_side2side);
  }
  NORETURN(Reals());
}

Reals get_edge_tangents(Mesh* mesh, LOs curv_edge2edge) {
  CHECK(mesh->dim() == 3);
  auto ncurv_edges = curv_edge2edge.size();
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto coords = mesh->coords();
  Write<Real> normals(ncurv_edges * 3);
  auto lambda = LAMBDA(LO curv_edge) {
    auto e = curv_edge2edge[curv_edge];
    auto v = gather_verts<2>(ev2v, e);
    auto x = gather_vectors<2, 3>(coords, v);
    auto b = simplex_basis<3, 1>(x);
    auto n = normalize(b[0]);
    set_vector(normals, curv_edge, n);
  };
  parallel_for(ncurv_edges, lambda);
  return normals;
}

Reals get_hinge_angles(Mesh* mesh, Reals surf_side_normals,
                       LOs surf_hinge2hinge, LOs side2surf_side) {
  switch (mesh->dim()) {
    case 3:
      return get_hinge_angles_tmpl<3>(mesh, surf_side_normals, surf_hinge2hinge,
                                      side2surf_side);
    case 2:
      return get_hinge_angles_tmpl<2>(mesh, surf_side_normals, surf_hinge2hinge,
                                      side2surf_side);
  }
  NORETURN(Reals());
}

/*
 * Max, Nelson.
 * "Weights for computing vertex normals from facet normals."
 * Journal of Graphics Tools 4.2 (1999): 1-6.
 */

static Reals tri_vert_normal_weights(Mesh* mesh, LOs surf_tri2tri) {
  auto ntris = mesh->ntris();
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

static Reals get_recip_length_weights(Mesh* mesh) {
  auto edge2len = measure_edges_real(mesh);
  auto nedges = mesh->nedges();
  auto edge2weight_w = Write<Real>(nedges);
  auto f = LAMBDA(LO se) { edge2weight_w[se] = 1.0 / edge2len[se]; };
  parallel_for(nedges, f);
  auto edge2weight = Reals(edge2weight_w);
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto ve2e = v2e.ab2b;
  auto ve2weight = unmap(ve2e, edge2weight, 1);
  return ve2weight;
}

static Reals side_vert_normal_weights(Mesh* mesh, LOs surf_side2side) {
  if (mesh->dim() == 3) return tri_vert_normal_weights(mesh, surf_side2side);
  if (mesh->dim() == 2) return get_recip_length_weights(mesh);
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

Reals get_vert_tangents(Mesh* mesh, LOs curv_edge2edge,
                        Reals curv_edge_tangents, LOs curv_vert2vert) {
  CHECK(mesh->dim() == 3);
  CHECK(mesh->owners_have_all_upward(VERT));
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto nve = ve2e.size();
  /* We would like to handle meshes where the mesh edges along
   * a model edge do not all point in the same direction
   * (this will be very common when we derive edges from only
   *  element connectivity).
   * In order to handle this case, we locally negate vectors as
   * necessary to get a correct tangent vector at each vertex
   */
  auto edge2curv_edge = invert_injective_map(curv_edge2edge, mesh->nedges());
  auto graph_tangents_w = Write<Real>(nve * 3);
  auto nv = mesh->nverts();
  auto neg_func = LAMBDA(LO v) {
    Int lc = 0;
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      auto e = ve2e[ve];
      auto ce = edge2curv_edge[e];
      if (-1 == ce) {
        set_vector(graph_tangents_w, ve, zero_vector<3>());
        continue;
      }
      auto code = v2e.codes[ve];
      auto eev = code_which_down(code);
      auto tangent = get_vector<3>(curv_edge_tangents, ce);
      if (eev == lc)
        set_vector(graph_tangents_w, ve, tangent);
      else
        set_vector(graph_tangents_w, ve, -tangent);
      ++lc;
    }
  };
  parallel_for(nv, neg_func);
  auto graph_tangents = Reals(graph_tangents_w);
  auto weights = get_recip_length_weights(mesh);
  auto vert_tangents =
      graph_weighted_average_arc_data(v2e, weights, graph_tangents, 3);
  vert_tangents = mesh->sync_array(VERT, vert_tangents, 3);
  auto curv_vert_tangents = unmap(curv_vert2vert, vert_tangents, 3);
  for (Int i = 0; i < curv_vert_tangents.size() / 3; ++i) {
    auto vec = get_vector<3>(curv_vert_tangents, i);
    fprintf(stderr, "vec[%d] = %f %f %f\n", i, vec[0], vec[1], vec[2]);
  }
  return normalize_vectors(curv_vert_tangents, 3);
}

}  // end namespace surf

}  // end namespace osh
