#include "Omega_h_confined.hpp"

#include "Omega_h_for.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_simplex.hpp"

/* Code to find classification-enforced constraints on
 * element size and shape.
 * For example, two geometric boundaries may be so close
 * as to force a mesh edge to be much smaller than desired,
 * similarly a sharp angle may be so acute that elements
 * contained by it are forced to be low-quality
 */

namespace Omega_h {

Bytes find_bridge_edges(Mesh* mesh) {
  auto verts2class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2class_dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto edges_are_bridges_w = Write<I8>(mesh->nedges());
  auto edges2verts = mesh->ask_verts_of(EDGE);
  auto f = OMEGA_H_LAMBDA(LO edge) {
    auto eev2v = gather_verts<2>(edges2verts, edge);
    auto edim = edges2class_dim[edge];
    edges_are_bridges_w[edge] = ((edim != verts2class_dim[eev2v[0]]) &&
                                 (edim != verts2class_dim[eev2v[1]]));
  };
  parallel_for(mesh->nedges(), f, "find_bridge_edges");
  return edges_are_bridges_w;
}

template <Int dim>
Reals get_edge_pad_dists(Mesh* mesh, Read<I8> edges_are_bridges) {
  auto coords = mesh->coords();
  auto edges2verts = mesh->ask_verts_of(EDGE);
  auto out = Write<Real>(mesh->nedges(), -1.0);
  auto f = OMEGA_H_LAMBDA(LO edge) {
    if (!edges_are_bridges[edge]) return;
    auto eev2v = gather_verts<2>(edges2verts, edge);
    auto eev2x = gather_vectors<2, dim>(coords, eev2v);
    auto h = norm(eev2x[1] - eev2x[0]);
    out[edge] = h;
  };
  parallel_for(mesh->nedges(), f);
  return out;
}

template <Int dim>
Reals get_tri_pad_dists(Mesh* mesh, Read<I8> edges_are_bridges) {
  auto coords = mesh->coords();
  auto tris2verts = mesh->ask_verts_of(FACE);
  auto tris2edges = mesh->ask_down(FACE, EDGE).ab2b;
  auto out = Write<Real>(mesh->nfaces(), -1.0);
  auto f = OMEGA_H_LAMBDA(LO tri) {
    auto ttv2v = gather_verts<3>(tris2verts, tri);
    auto ttv2x = gather_vectors<3, dim>(coords, ttv2v);
    auto tte2e = gather_down<3>(tris2edges, tri);
    auto tte2b = gather_values<3>(edges_are_bridges, tte2e);
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
      out[tri] = h;
    }
  };
  parallel_for(mesh->nfaces(), f);
  return out;
}

Reals get_tet_pad_dists(Mesh* mesh, Read<I8> edges_are_bridges) {
  auto coords = mesh->coords();
  auto tets2verts = mesh->ask_verts_of(REGION);
  auto tets2edges = mesh->ask_down(REGION, EDGE).ab2b;
  auto out = Write<Real>(mesh->nregions(), -1.0);
  auto f = OMEGA_H_LAMBDA(LO tet) {
    auto ttv2v = gather_verts<4>(tets2verts, tet);
    auto ttv2x = gather_vectors<4, 3>(coords, ttv2v);
    auto tte2e = gather_down<6>(tets2edges, tet);
    auto tte2b = gather_values<6>(edges_are_bridges, tte2e);
    auto nb = reduce(tte2b, plus<I8>());
    if (nb == 0) return;
    if (nb == 4) {
      for (Int tte = 0; tte < 3; ++tte) {
        if (tte2b[tte]) continue;
        auto opp = simplex_opposite_template(REGION, EDGE, tte);
        if (tte2b[opp]) continue;
        // at this point we have edge-edge nearness
        auto a = ttv2x[simplex_down_template(REGION, EDGE, tte, 0)];
        auto b = ttv2x[simplex_down_template(REGION, EDGE, tte, 1)];
        auto c = ttv2x[simplex_down_template(REGION, EDGE, opp, 0)];
        auto d = ttv2x[simplex_down_template(REGION, EDGE, opp, 1)];
        auto ab = b - a;
        auto cd = d - c;
        auto n = normalize(cross(ab, cd));
        auto h = (a - c) * n;
        // project onto the normal plane
        a = a - (n * (n * a));
        b = b - (n * (n * b));
        c = c - (n * (n * c));
        d = d - (n * (n * d));
        if (!((get_triangle_normal(a, b, c) * get_triangle_normal(a, b, d)) <
                    0 &&
                (get_triangle_normal(c, d, a) * get_triangle_normal(c, d, b)) <
                    0)) {
          break;
        }
        out[tet] = h;
        return;  // edge-edge implies no plane-vertex
      }
    }
    Real h = ArithTraits<Real>::max();  // multiple vertex-planes may occur
    for (Int ttv = 0; ttv < 4; ++ttv) {
      Few<Int, 3> vve2tte;
      Few<Int, 3> vve2wd;
      for (Int vve = 0; vve < 3; ++vve) {
        vve2tte[vve] = simplex_up_template(REGION, VERT, ttv, vve).up;
        vve2wd[vve] = simplex_up_template(REGION, VERT, ttv, vve).which_down;
      }
      Few<I8, 3> vve2b;
      for (Int vve = 0; vve < 3; ++vve) vve2b[vve] = tte2b[vve2tte[vve]];
      if (reduce(vve2b, plus<I8>()) != 3) continue;
      // at this point, we have vertex-plane nearness
      auto o = ttv2x[ttv];
      Few<Int, 3> vve2ttv;
      for (Int vve = 0; vve < 3; ++vve) {
        vve2ttv[vve] =
            simplex_down_template(REGION, EDGE, vve2tte[vve], 1 - vve2wd[vve]);
      }
      Few<Vector<3>, 3> vve2x;
      for (Int vve = 0; vve < 3; ++vve) vve2x[vve] = ttv2x[vve2ttv[vve]];
      auto a = vve2x[0];
      auto b = vve2x[1];
      auto c = vve2x[2];
      auto ab = b - a;
      auto ac = c - a;
      auto n = normalize(cross(ab, ac));
      auto oa = a - o;
      auto od = n * (n * oa);
      Matrix<3, 2> basis;
      basis[0] = ab;
      basis[1] = ac;
      auto inv_basis = pseudo_invert(basis);
      auto ad = od - oa;
      auto xi = form_barycentric(inv_basis * ad);
      if (!is_barycentric_inside(xi)) continue;
      auto this_h = norm(od);
      h = min2(h, this_h);
    }
    if (h == ArithTraits<Real>::max()) h = -1.0;
    out[tet] = h;
  };
  parallel_for(mesh->nregions(), f);
  return out;
}

Reals get_pad_dists(Mesh* mesh, Int pad_dim, Read<I8> edges_are_bridges) {
  if (pad_dim == EDGE) {
    if (mesh->dim() == 3) {
      return get_edge_pad_dists<3>(mesh, edges_are_bridges);
    } else if (mesh->dim() == 2) {
      return get_edge_pad_dists<2>(mesh, edges_are_bridges);
    } else if (mesh->dim() == 1) {
      return get_edge_pad_dists<1>(mesh, edges_are_bridges);
    }
  } else if (pad_dim == FACE) {
    if (mesh->dim() == 3) {
      return get_tri_pad_dists<3>(mesh, edges_are_bridges);
    } else if (mesh->dim() == 2) {
      return get_tri_pad_dists<2>(mesh, edges_are_bridges);
    }
  } else if (pad_dim == REGION) {
    return get_tet_pad_dists(mesh, edges_are_bridges);
  }
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals get_pinched_tri_angles_dim(Mesh* mesh) {
  auto verts2class_dim = mesh->get_array<I8>(VERT, "class_dim");
  auto edges2class_dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto tris2edges = mesh->ask_down(FACE, EDGE).ab2b;
  auto tris2verts = mesh->ask_down(FACE, VERT).ab2b;
  auto coords = mesh->coords();
  auto tri_angles_w = Write<Real>(mesh->nfaces());
  auto f = OMEGA_H_LAMBDA(LO tri) {
    auto ttv2v = gather_down<3>(tris2verts, tri);
    auto tte2e = gather_down<3>(tris2edges, tri);
    auto ttv2dim = gather_values(verts2class_dim, ttv2v);
    auto tte2dim = gather_values(edges2class_dim, tte2e);
    auto ttv2x = gather_vectors<3, dim>(coords, ttv2v);
    Real tri_angle = -1.0;
    for (Int i = 0; i < 3; ++i) {
      if (ttv2dim[i] != 0) continue;
      if (tte2dim[(i + 0) % 3] != 1) continue;
      if (tte2dim[(i + 2) % 3] != 1) continue;
      auto d0 = ttv2x[(i + 1) % 3] - ttv2x[i];
      auto d1 = ttv2x[(i + 2) % 3] - ttv2x[i];
      auto cos_a = normalize(d0) * normalize(d1);
      auto angle = std::acos(cos_a);
      if (tri_angle == -1.0 || angle < tri_angle) {
        tri_angle = angle;
      }
    }
    tri_angles_w[tri] = tri_angle;
  };
  parallel_for(mesh->nfaces(), f, "get_pinched_tri_angles");
  return tri_angles_w;
}

Reals get_pinched_tet_angles(Mesh* mesh) {
  auto edges2class_dim = mesh->get_array<I8>(EDGE, "class_dim");
  auto tris2class_dim = mesh->get_array<I8>(FACE, "class_dim");
  auto tets2edges = mesh->ask_down(REGION, EDGE).ab2b;
  auto tets2tris = mesh->ask_down(REGION, FACE).ab2b;
  auto edges2verts = mesh->ask_down(EDGE, VERT).ab2b;
  auto coords = mesh->coords();
  auto tet_angles_w = Write<Real>(mesh->nregions());
  auto f = OMEGA_H_LAMBDA(LO tet) {
    auto kke2e = gather_down<6>(tets2edges, tet);
    auto kkt2t = gather_down<4>(tets2tris, tet);
    auto kke2dim = gather_values(edges2class_dim, kke2e);
    auto kkt2dim = gather_values(tris2class_dim, kkt2t);
    Real tet_angle = -1.0;
    for (Int kke = 0; kke < 6; ++kke) {
      if (kke2dim[kke] != 1) continue;
      auto ket0 = simplex_up_template(3, 1, kke, 0).up;
      auto ket1 = simplex_up_template(3, 1, kke, 1).up;
      if (kkt2dim[ket0] != 2) continue;
      if (kkt2dim[ket1] != 2) continue;
      auto e0 = kke2e[kke];
      auto e0v2v = gather_down<2>(edges2verts, e0);
      auto e0v2x = gather_vectors<2, 3>(coords, e0v2v);
      auto tangent = e0v2x[1] - e0v2x[0];
      auto kke_opp = simplex_opposite_template(3, 1, kke);
      auto e1 = kke2e[kke_opp];
      auto e1v2v = gather_down<2>(edges2verts, e1);
      auto e1v2x = gather_vectors<2, 3>(coords, e1v2v);
      auto l0 = e1v2x[0] - e0v2x[0];
      auto l1 = e1v2x[1] - e0v2x[0];
      auto d0 = reject(l0, tangent);
      auto d1 = reject(l1, tangent);
      auto cos_a = normalize(d0) * normalize(d1);
      auto angle = std::acos(cos_a);
      if (tet_angle == -1.0 || angle < tet_angle) {
        tet_angle = angle;
      }
    }
    tet_angles_w[tet] = tet_angle;
  };
  parallel_for(mesh->nregions(), f, "get_pinched_tet_angles");
  return tet_angles_w;
}

Reals get_pinched_angles(Mesh* mesh, Int pinched_dim) {
  if (pinched_dim == FACE) {
    if (mesh->dim() == 3) {
      return get_pinched_tri_angles_dim<3>(mesh);
    } else if (mesh->dim() == 2) {
      return get_pinched_tri_angles_dim<2>(mesh);
    }
  } else if (pinched_dim == REGION) {
    return get_pinched_tet_angles(mesh);
  }
  OMEGA_H_NORETURN(Reals());
}

}  // end namespace Omega_h
