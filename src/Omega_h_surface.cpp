#include "Omega_h_surface.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_shape.hpp"

namespace Omega_h {

namespace {

Reals get_triangle_normals(Mesh* mesh, LOs surf_tri2tri) {
  OMEGA_H_CHECK(mesh->dim() == 3);
  auto nsurf_tris = surf_tri2tri.size();
  auto fv2v = mesh->ask_verts_of(FACE);
  auto coords = mesh->coords();
  Write<Real> normals(nsurf_tris * 3);
  auto lambda = OMEGA_H_LAMBDA(LO surf_tri) {
    auto f = surf_tri2tri[surf_tri];
    auto v = gather_verts<3>(fv2v, f);
    auto x = gather_vectors<3, 3>(coords, v);
    auto b = simplex_basis<3, 2>(x);
    auto n = normalize(cross(b[0], b[1]));
    set_vector(normals, surf_tri, n);
  };
  parallel_for(nsurf_tris, lambda, "get_triangle_normals");
  return normals;
}

Reals get_edge_normals(Mesh* mesh, LOs surf_edge2edge) {
  OMEGA_H_CHECK(mesh->dim() == 2);
  auto nsurf_edges = surf_edge2edge.size();
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto coords = mesh->coords();
  Write<Real> normals(nsurf_edges * 2);
  auto lambda = OMEGA_H_LAMBDA(LO surf_edge) {
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
  parallel_for(nsurf_edges, lambda, "get_edge_normals");
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
  auto f = OMEGA_H_LAMBDA(LO surf_hinge) {
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
    angles[surf_hinge] = std::acos(n[0] * n[1]);
  };
  parallel_for(nsurf_hinges, f, "get_hingle_angles");
  return angles;
}

}  // end anonymous namespace

Reals get_side_vectors(Mesh* mesh, LOs surf_side2side) {
  if (mesh->dim() == 3) return get_triangle_normals(mesh, surf_side2side);
  if (mesh->dim() == 2) return get_edge_normals(mesh, surf_side2side);
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals get_curv_edge_tangents_dim(Mesh* mesh, LOs curv_edge2edge) {
  OMEGA_H_CHECK(mesh->dim() == dim);
  auto ncurv_edges = curv_edge2edge.size();
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto coords = mesh->coords();
  Write<Real> normals(ncurv_edges * dim);
  auto lambda = OMEGA_H_LAMBDA(LO curv_edge) {
    auto e = curv_edge2edge[curv_edge];
    auto v = gather_verts<2>(ev2v, e);
    auto x = gather_vectors<2, dim>(coords, v);
    auto b = simplex_basis<dim, 1>(x);
    auto n = normalize(b[0]);
    set_vector(normals, curv_edge, n);
  };
  parallel_for(ncurv_edges, lambda, "get_curv_edge_tangents");
  return normals;
}

Reals get_curv_edge_tangents(Mesh* mesh, LOs curv_edge2edge) {
  if (mesh->dim() == 3)
    return get_curv_edge_tangents_dim<3>(mesh, curv_edge2edge);
  if (mesh->dim() == 2)
    return get_curv_edge_tangents_dim<2>(mesh, curv_edge2edge);
  OMEGA_H_NORETURN(Reals());
}

Reals get_hinge_angles(Mesh* mesh, Reals surf_side_normals,
    LOs surf_hinge2hinge, LOs side2surf_side) {
  switch (mesh->dim()) {
    case 3:
      return get_hinge_angles_tmpl<3>(
          mesh, surf_side_normals, surf_hinge2hinge, side2surf_side);
    case 2:
      return get_hinge_angles_tmpl<2>(
          mesh, surf_side_normals, surf_hinge2hinge, side2surf_side);
  }
  OMEGA_H_NORETURN(Reals());
}

/*
 * Max, Nelson.
 * "Weights for computing vertex normals from facet normals."
 * Journal of Graphics Tools 4.2 (1999): 1-6.
 */

Reals tri_vert_normal_weights(
    Mesh* mesh, Adj surf_verts2tris, LOs tri2surf_tri) {
  auto nsurf_verts = surf_verts2tris.nnodes();
  auto fv2v = mesh->ask_verts_of(FACE);
  auto coords = mesh->coords();
  auto nweights = surf_verts2tris.nedges();
  auto weights = Write<Real>(nweights);
  auto func = OMEGA_H_LAMBDA(LO surf_vert) {
    for (auto vf = surf_verts2tris.a2ab[surf_vert];
         vf < surf_verts2tris.a2ab[surf_vert + 1]; ++vf) {
      auto f = surf_verts2tris.ab2b[vf];
      if (tri2surf_tri[f] < 0) {
        weights[vf] = 0.0;
        continue;
      }
      auto code = surf_verts2tris.codes[vf];
      auto ffv = code_which_down(code);
      auto rot = rotation_to_first(3, ffv);
      auto ffv2v_0 = gather_verts<3>(fv2v, f);
      Few<LO, 3> ffv2v;
      rotate_adj<3>(rot, ffv2v_0, 0, ffv2v, 0);
      auto ffv2p = gather_vectors<3, 3>(coords, ffv2v);
      auto b = simplex_basis<3, 2>(ffv2p);
      auto w =
          norm(cross(b[0], b[1])) / (norm_squared(b[0]) * norm_squared(b[1]));
      weights[vf] = w;
    }
  };
  parallel_for(nsurf_verts, func, "tri_vert_normal_weights");
  return weights;
}

/* following the same derivation as Max did above,
 * we find that the weights for edges should be
 * one over the length of the edge.
 */

Reals get_recip_length_weights(
    Mesh* mesh, Adj surf_verts2edges, LOs edge2surf_edge) {
  auto edge2len = measure_edges_real(mesh);
  auto nedges = mesh->nedges();
  auto edge2weight_w = Write<Real>(nedges);
  auto f = OMEGA_H_LAMBDA(LO edge) {
    if (edge2surf_edge[edge] == -1)
      edge2weight_w[edge] = 0.0;
    else
      edge2weight_w[edge] = 1.0 / edge2len[edge];
  };
  parallel_for(nedges, f, "get_recip_length_weights");
  auto edge2weight = Reals(edge2weight_w);
  return unmap(surf_verts2edges.ab2b, edge2weight, 1);
}

static Reals side_vert_normal_weights(
    Mesh* mesh, Adj surf_verts2sides, LOs side2surf_side) {
  if (mesh->dim() == 3) {
    return tri_vert_normal_weights(mesh, surf_verts2sides, side2surf_side);
  }
  if (mesh->dim() == 2) {
    return get_recip_length_weights(mesh, surf_verts2sides, side2surf_side);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals get_side_vert_normals(Mesh* mesh, LOs surf_side2side,
    Reals surf_side_normals, LOs surf_vert2vert) {
  OMEGA_H_CHECK(mesh->owners_have_all_upward(VERT));
  auto dim = mesh->dim();
  auto verts2sides = mesh->ask_up(VERT, dim - 1);
  auto surf_verts2sides = unmap_adjacency(surf_vert2vert, verts2sides);
  auto nsides = mesh->nents(dim - 1);
  auto side2surf_side = invert_injective_map(surf_side2side, nsides);
  auto weights =
      side_vert_normal_weights(mesh, surf_verts2sides, side2surf_side);
  auto side_normals =
      map_onto(surf_side_normals, surf_side2side, nsides, 0.0, dim);
  auto surf_vert_normals =
      graph_weighted_average(surf_verts2sides, weights, side_normals, dim);
  surf_vert_normals = mesh->sync_subset_array(
      VERT, surf_vert_normals, surf_vert2vert, 0.0, dim);
  return normalize_vectors(surf_vert_normals, dim);
}

/* We would like to handle meshes where the mesh edges along
 * a model edge do not all point in the same direction
 * (this will be very common when we derive edges from only
 *  element connectivity).
 * In order to handle this case, we locally negate vectors as
 * necessary to get a correct tangent vector at each vertex
 */
Read<I8> get_curv_vert_edge_flips(
    Adj curv_verts2edges, LOs edge2curv_edge) {
  auto out = Write<I8>(curv_verts2edges.nedges());
  auto ncurv_verts = curv_verts2edges.nnodes();
  auto f = OMEGA_H_LAMBDA(LO curv_vert) {
    Int lc = 0;
    for (auto ve = curv_verts2edges.a2ab[curv_vert];
         ve < curv_verts2edges.a2ab[curv_vert + 1]; ++ve) {
      auto e = curv_verts2edges.ab2b[ve];
      auto ce = edge2curv_edge[e];
      if (-1 == ce) continue;
      auto code = curv_verts2edges.codes[ve];
      auto eev = code_which_down(code);
      out[ve] = I8(eev != lc);
      ++lc;
    }
  };
  parallel_for(ncurv_verts, f, "get_curv_vert_edge_flips");
  return out;
}

Read<I8> get_curv_edge_vert_flips(
    Mesh* mesh, Adj curv_verts2edges, LOs edge2curv_edge) {
  auto in = get_curv_vert_edge_flips(curv_verts2edges, edge2curv_edge);
  auto out = Write<I8>(mesh->nedges() * 2, I8(0));
  auto ncurv_verts = curv_verts2edges.nnodes();
  auto f = OMEGA_H_LAMBDA(LO curv_vert) {
    for (auto ve = curv_verts2edges.a2ab[curv_vert];
         ve < curv_verts2edges.a2ab[curv_vert + 1]; ++ve) {
      auto e = curv_verts2edges.ab2b[ve];
      auto code = curv_verts2edges.codes[ve];
      auto eev = code_which_down(code);
      out[e * 2 + eev] = in[ve];
    }
  };
  parallel_for(ncurv_verts, f, "get_curv_edge_vert_flips");
  return out;
}

template <Int dim>
Reals get_curv_vert_tangents_dim(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_tangents, LOs curv_vert2vert) {
  OMEGA_H_CHECK(mesh->dim() == dim);
  OMEGA_H_CHECK(mesh->owners_have_all_upward(VERT));
  auto verts2edges = mesh->ask_up(VERT, EDGE);
  auto curv_verts2edges = unmap_adjacency(curv_vert2vert, verts2edges);
  auto narcs = curv_verts2edges.nedges();
  auto edge2curv_edge = invert_injective_map(curv_edge2edge, mesh->nedges());
  auto arcs2tangents_w = Write<Real>(narcs * dim);
  auto arcs2flip = get_curv_vert_edge_flips(curv_verts2edges, edge2curv_edge);
  auto f = OMEGA_H_LAMBDA(LO arc) {
    auto e = curv_verts2edges.ab2b[arc];
    auto ce = edge2curv_edge[e];
    if (-1 == ce)
      set_vector(arcs2tangents_w, arc, zero_vector<dim>());
    else {
      auto tangent = get_vector<dim>(curv_edge_tangents, ce);
      if (arcs2flip[arc]) tangent = -tangent;
      set_vector(arcs2tangents_w, arc, tangent);
    }
  };
  parallel_for(narcs, f, "get_curv_vert_tangents");
  auto arcs2tangents = Reals(arcs2tangents_w);
  auto weights =
      get_recip_length_weights(mesh, curv_verts2edges, edge2curv_edge);
  auto curv_vert_tangents = graph_weighted_average_arc_data(
      curv_verts2edges, weights, arcs2tangents, dim);
  curv_vert_tangents = mesh->sync_subset_array(
      VERT, curv_vert_tangents, curv_vert2vert, 0.0, dim);
  return normalize_vectors(curv_vert_tangents, dim);
}

Reals get_curv_vert_tangents(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_tangents, LOs curv_vert2vert) {
  if (mesh->dim() == 3) {
    return get_curv_vert_tangents_dim<3>(
        mesh, curv_edge2edge, curv_edge_tangents, curv_vert2vert);
  }
  if (mesh->dim() == 2) {
    return get_curv_vert_tangents_dim<2>(
        mesh, curv_edge2edge, curv_edge_tangents, curv_vert2vert);
  }
  OMEGA_H_NORETURN(Reals());
}

/* Rusinkiewicz, Szymon.
 * "Estimating curvatures and their derivatives on triangle meshes."
 * 3D Data Processing, Visualization and Transmission, 2004.
 * 3DPVT 2004. Proceedings. 2nd International Symposium on. IEEE, 2004.
 */

Reals get_surf_tri_IIs(Mesh* mesh, LOs surf_tri2tri, Reals surf_tri_normals,
    LOs surf_vert2vert, Reals surf_vert_normals) {
  auto vert2surf_vert = invert_injective_map(surf_vert2vert, mesh->nverts());
  auto tris2verts = mesh->ask_verts_of(FACE);
  auto coords = mesh->coords();
  auto nsurf_tris = surf_tri2tri.size();
  auto surf_tri_IIs_w = Write<Real>(nsurf_tris * symm_ncomps(2));
  auto f = OMEGA_H_LAMBDA(LO surf_tri) {
    auto tri = surf_tri2tri[surf_tri];
    auto tn = get_vector<3>(surf_tri_normals, surf_tri);
    auto nuv = form_ortho_basis(tn);
    auto u = nuv[1];
    auto v = nuv[2];
    auto ttv = gather_verts<3>(tris2verts, tri);
    auto p = gather_vectors<3, 3>(coords, ttv);
    Few<Vector<3>, 3> n;
    for (Int i = 0; i < 3; ++i) {
      auto vert = ttv[i];
      auto surf_vert = vert2surf_vert[vert];
      if (surf_vert >= 0)
        n[i] = get_vector<3>(surf_vert_normals, surf_vert);
      else
        n[i] = tn;
    }
    Few<Vector<3>, 3> e;
    for (Int i = 0; i < 3; ++i) e[i] = p[(i + 2) % 3] - p[(i + 1) % 3];
    Few<Vector<3>, 3> dn;
    for (Int i = 0; i < 3; ++i) dn[i] = n[(i + 2) % 3] - n[(i + 1) % 3];
    Matrix<6, 3> A;
    Vector<6> rhs;
    for (Int i = 0; i < 3; ++i) {
      A[2][i * 2 + 1] = A[0][i * 2 + 0] = e[i] * u;
      A[1][i * 2 + 1] = A[2][i * 2 + 0] = e[i] * v;
      A[1][i * 2 + 0] = A[0][i * 2 + 1] = 0;
      rhs[i * 2 + 0] = dn[i] * u;
      rhs[i * 2 + 1] = dn[i] * v;
    }
    auto II_comps = solve_using_qr(A, rhs);
    auto II = vector2symm(II_comps);
    set_symm(surf_tri_IIs_w, surf_tri, II);
  };
  parallel_for(nsurf_tris, f, "get_surf_tri_IIs");
  return surf_tri_IIs_w;
}

/* Meyer, Mark, et al.
 * "Discrete differential-geometry operators for triangulated 2-manifolds."
 * Visualization and mathematics III. Springer Berlin Heidelberg, 2003. 35-57.
 */

OMEGA_H_INLINE Real get_mixed_area(Few<Vector<3>, 3> p, Int ttv) {
  Few<Vector<3>, 3> e;
  for (Int i = 0; i < 3; ++i) e[i] = p[(i + 1) % 3] - p[i];
  Int ttv_obtuse = -1;
  for (Int i = 0; i < 3; ++i) {
    if (e[i] * e[(i + 2) % 3] > 0.0) {
      ttv_obtuse = i;
      break;
    }
  }
  Few<Vector<3>, 2> basis;
  basis[0] = e[ttv];
  basis[1] = -e[(ttv + 2) % 3];
  if (ttv_obtuse >= 0) {
    auto area = triangle_area_from_basis(basis);
    if (ttv == ttv_obtuse)
      return area / 2.0;
    else
      return area / 4.0;
  } else {
    auto ao = get_circumcenter_vector(basis);
    Few<Vector<3>, 2> halfbasis;
    halfbasis[0] = basis[0] / 2.0;
    halfbasis[1] = ao;
    auto area = triangle_area_from_basis(halfbasis);
    halfbasis[0] = ao;
    halfbasis[1] = basis[1] / 2.0;
    area += triangle_area_from_basis(halfbasis);
    return area;
  }
}

OMEGA_H_INLINE Matrix<3, 3> rotate_to_plane(Vector<3> n, Matrix<3, 3> tnuv) {
  auto tn = tnuv[0];
  auto cp = cross(tn, n);
  auto cpl = norm(cp);
  if (cpl < EPSILON) return tnuv;
  auto axis = cp / cpl;
  auto angle = std::asin(cpl);
  auto r = rotate(angle, axis);
  return r * tnuv;
}

Reals get_surf_vert_IIs(Mesh* mesh, LOs surf_tri2tri, Reals surf_tri_normals,
    Reals surf_tri_IIs, LOs surf_vert2vert, Reals surf_vert_normals) {
  auto nsurf_verts = surf_vert2vert.size();
  auto verts2tris = mesh->ask_up(VERT, FACE);
  auto tri2surf_tri = invert_injective_map(surf_tri2tri, mesh->nfaces());
  auto coords = mesh->coords();
  auto tris2verts = mesh->ask_verts_of(FACE);
  auto verts_not_surf = map_onto(
      Read<I8>(nsurf_verts, I8(0)), surf_vert2vert, mesh->nverts(), I8(1), 1);
  auto tris_touch_bdry = mark_up(mesh, VERT, FACE, verts_not_surf);
  auto surf_vert_IIs_w = Write<Real>(nsurf_verts * 3);
  auto f = OMEGA_H_LAMBDA(LO surf_vert) {
    auto vert = surf_vert2vert[surf_vert];
    auto n = get_vector<3>(surf_vert_normals, surf_vert);
    auto nuv = form_ortho_basis(n);
    Real ws = 0.0;
    Vector<3> comps = vector_3(0.0, 0.0, 0.0);
    Int nadj_int_tris = 0;
    for (auto vt = verts2tris.a2ab[vert]; vt < verts2tris.a2ab[vert + 1];
         ++vt) {
      auto tri = verts2tris.ab2b[vt];
      auto surf_tri = tri2surf_tri[tri];
      if (surf_tri < 0) continue;
      if (!tris_touch_bdry[tri]) ++nadj_int_tris;
    }
    for (auto vt = verts2tris.a2ab[vert]; vt < verts2tris.a2ab[vert + 1];
         ++vt) {
      auto tri = verts2tris.ab2b[vt];
      auto surf_tri = tri2surf_tri[tri];
      if (surf_tri < 0) continue;
      if (nadj_int_tris && tris_touch_bdry[tri]) continue;
      auto tn = get_vector<3>(surf_tri_normals, surf_tri);
      auto fnuv = form_ortho_basis(tn);
      fnuv = rotate_to_plane(n, fnuv);
      Matrix<2, 2> jac;
      for (Int i = 0; i < 2; ++i) {
        for (Int j = 0; j < 2; ++j) {
          jac[j][i] = nuv[j + 1] * fnuv[i + 1];
        }
      }
      auto tri_II = get_symm<2>(surf_tri_IIs, surf_tri);
      Vector<3> tri_comps;
      tri_comps[0] = jac[0] * (tri_II * jac[0]);
      tri_comps[1] = jac[1] * (tri_II * jac[1]);
      tri_comps[2] = jac[0] * (tri_II * jac[1]);
      auto code = verts2tris.codes[vt];
      auto ttv = code_which_down(code);
      auto ttv2v = gather_verts<3>(tris2verts, tri);
      auto p = gather_vectors<3, 3>(coords, ttv2v);
      auto w = get_mixed_area(p, ttv);
      comps = comps + (tri_comps * w);
      ws += w;
    }
    comps = comps / ws;
    set_vector(surf_vert_IIs_w, surf_vert, comps);
  };
  parallel_for(nsurf_verts, f, "get_surf_vert_IIs");
  auto surf_vert_IIs = Reals(surf_vert_IIs_w);
  return mesh->sync_subset_array(VERT, surf_vert_IIs, surf_vert2vert, 0.0, 3);
}

template <Int dim>
Reals get_curv_edge_curvatures_dim(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_tangents, LOs curv_vert2vert, Reals curv_vert_tangents) {
  auto vert2curv_vert = invert_injective_map(curv_vert2vert, mesh->nverts());
  auto edge2curv_edge = invert_injective_map(curv_edge2edge, mesh->nedges());
  auto edges2verts = mesh->ask_verts_of(EDGE);
  auto ncurv_edges = curv_edge2edge.size();
  auto curv_edges2curvature_w = Write<Real>(ncurv_edges);
  auto verts2edges = mesh->ask_up(VERT, EDGE);
  auto curv_verts2edges = unmap_adjacency(curv_vert2vert, verts2edges);
  auto ev_flips =
      get_curv_edge_vert_flips(mesh, curv_verts2edges, edge2curv_edge);
  auto coords = mesh->coords();
  auto f = OMEGA_H_LAMBDA(LO curv_edge) {
    auto edge = curv_edge2edge[curv_edge];
    auto et = get_vector<dim>(curv_edge_tangents, curv_edge);
    Few<Vector<dim>, 2> ts;
    auto eev2v = gather_verts<2>(edges2verts, edge);
    for (Int eev = 0; eev < 2; ++eev) {
      auto vert = eev2v[eev];
      auto curv_vert = vert2curv_vert[vert];
      Vector<dim> vt;
      if (curv_vert >= 0) {
        vt = get_vector<dim>(curv_vert_tangents, curv_vert);
        if (ev_flips[edge * 2 + eev]) vt = -vt;
      } else {
        vt = et;
      }
      ts[eev] = vt;
    }
    auto p = gather_vectors<2, dim>(coords, eev2v);
    auto e = p[1] - p[0];
    auto l = norm(e);
    auto u = e / l;
    auto dt = ts[1] - ts[0];
    auto curvature = norm(dt - (u * (dt * u))) / l;
    curv_edges2curvature_w[curv_edge] = curvature;
  };
  parallel_for(ncurv_edges, f, "get_curv_edge_curvatures");
  auto curv_edges2curvature = Reals(curv_edges2curvature_w);
  return mesh->sync_subset_array(
      EDGE, curv_edges2curvature, curv_edge2edge, 0.0, 1);
}

Reals get_curv_edge_curvatures(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_tangents, LOs curv_vert2vert, Reals curv_vert_tangents) {
  if (mesh->dim() == 3) {
    return get_curv_edge_curvatures_dim<3>(mesh, curv_edge2edge,
        curv_edge_tangents, curv_vert2vert, curv_vert_tangents);
  }
  if (mesh->dim() == 2) {
    return get_curv_edge_curvatures_dim<2>(mesh, curv_edge2edge,
        curv_edge_tangents, curv_vert2vert, curv_vert_tangents);
  }
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals get_curv_vert_curvatures_dim(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_curvatures, LOs curv_vert2vert) {
  auto verts2edges = mesh->ask_up(VERT, EDGE);
  auto edges2verts = mesh->ask_verts_of(EDGE);
  auto edge2curv_edge = invert_injective_map(curv_edge2edge, mesh->nedges());
  auto ncurv_verts = curv_vert2vert.size();
  auto verts_not_curv = map_onto(
      Read<I8>(ncurv_verts, I8(0)), curv_vert2vert, mesh->nverts(), I8(1), 1);
  auto edges_touch_bdry = mark_up(mesh, VERT, EDGE, verts_not_curv);
  auto coords = mesh->coords();
  auto curv_vert_curvatures_w = Write<Real>(ncurv_verts);
  auto f = OMEGA_H_LAMBDA(LO curv_vert) {
    auto vert = curv_vert2vert[curv_vert];
    Int nadj_int_edges = 0;
    for (auto ve = verts2edges.a2ab[vert]; ve < verts2edges.a2ab[vert + 1];
         ++ve) {
      auto edge = verts2edges.ab2b[ve];
      nadj_int_edges += !edges_touch_bdry[edge];
    }
    Real ws = 0.0;
    Real curvature = 0.0;
    for (auto ve = verts2edges.a2ab[vert]; ve < verts2edges.a2ab[vert + 1];
         ++ve) {
      auto edge = verts2edges.ab2b[ve];
      auto curv_edge = edge2curv_edge[edge];
      if (curv_edge < 0) continue;
      if (nadj_int_edges && edges_touch_bdry[edge]) continue;
      auto eev2v = gather_verts<2>(edges2verts, edge);
      auto p = gather_vectors<2, dim>(coords, eev2v);
      auto l = norm(p[1] - p[0]);
      auto ec = curv_edge_curvatures[curv_edge];
      curvature += ec * l;
      ws += l;
    }
    curvature /= ws;
    curv_vert_curvatures_w[curv_vert] = curvature;
  };
  parallel_for(ncurv_verts, f, "get_curv_vert_curvatures");
  auto curv_vert_curvatures = Reals(curv_vert_curvatures_w);
  return mesh->sync_subset_array(
      VERT, curv_vert_curvatures, curv_vert2vert, 0.0, 1);
}

Reals get_curv_vert_curvatures(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_curvatures, LOs curv_vert2vert) {
  if (mesh->dim() == 3) {
    return get_curv_vert_curvatures_dim<3>(
        mesh, curv_edge2edge, curv_edge_curvatures, curv_vert2vert);
  }
  if (mesh->dim() == 2) {
    return get_curv_vert_curvatures_dim<2>(
        mesh, curv_edge2edge, curv_edge_curvatures, curv_vert2vert);
  }
  OMEGA_H_NORETURN(Reals());
}

SurfaceInfo get_surface_info(Mesh* mesh) {
  SurfaceInfo out;
  if (mesh->dim() == 1) return out;
  auto sdim = mesh->dim() - 1;
  auto sides_are_surf = mark_by_class_dim(mesh, sdim, sdim);
  auto verts_are_surf = mark_by_class_dim(mesh, VERT, sdim);
  auto surf_side2side = collect_marked(sides_are_surf);
  auto surf_vert2vert = collect_marked(verts_are_surf);
  LOs curv_edge2edge;
  LOs curv_vert2vert;
  if (mesh->dim() == 3) {
    auto surf_side_normals = get_side_vectors(mesh, surf_side2side);
    auto surf_vert_normals = get_side_vert_normals(
        mesh, surf_side2side, surf_side_normals, surf_vert2vert);
    auto edges_are_curv = mark_by_class_dim(mesh, EDGE, EDGE);
    auto verts_are_curv = mark_by_class_dim(mesh, VERT, EDGE);
    curv_edge2edge = collect_marked(edges_are_curv);
    curv_vert2vert = collect_marked(verts_are_curv);
    auto surf_tri_IIs = get_surf_tri_IIs(mesh, surf_side2side,
        surf_side_normals, surf_vert2vert, surf_vert_normals);
    auto surf_vert_IIs = get_surf_vert_IIs(mesh, surf_side2side,
        surf_side_normals, surf_tri_IIs, surf_vert2vert, surf_vert_normals);
    out.surf_vert2vert = surf_vert2vert;
    out.surf_vert_normals = surf_vert_normals;
    out.surf_vert_IIs = surf_vert_IIs;
  } else {
    curv_edge2edge = surf_side2side;
    curv_vert2vert = surf_vert2vert;
  }
  auto curv_edge_tangents = get_curv_edge_tangents(mesh, curv_edge2edge);
  auto curv_vert_tangents = get_curv_vert_tangents(
      mesh, curv_edge2edge, curv_edge_tangents, curv_vert2vert);
  auto curv_edge_curvatures = get_curv_edge_curvatures(mesh, curv_edge2edge,
      curv_edge_tangents, curv_vert2vert, curv_vert_tangents);
  auto curv_vert_curvatures = get_curv_vert_curvatures(
      mesh, curv_edge2edge, curv_edge_curvatures, curv_vert2vert);
  out.curv_vert2vert = curv_vert2vert;
  out.curv_vert_tangents = curv_vert_tangents;
  out.curv_vert_curvatures = curv_vert_curvatures;
  return out;
}

Reals get_vert_curvatures(Mesh* mesh, SurfaceInfo surface_info) {
  Write<Real> out(mesh->nverts(), 0.0);
  if (mesh->dim() >= 3) {
    auto surf_vert_curvatures =
        get_max_eigenvalues(2, surface_info.surf_vert_IIs);
    map_into(surf_vert_curvatures, surface_info.surf_vert2vert, out, 1);
  }
  if (mesh->dim() >= 2) {
    map_into(
        surface_info.curv_vert_curvatures, surface_info.curv_vert2vert, out, 1);
  }
  return out;
}

}  // end namespace Omega_h
