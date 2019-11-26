#include "Omega_h_coarsen.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_collapse.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_most_normal.hpp"
#include "Omega_h_shape.hpp"

namespace Omega_h {

template <Int dim>
Reals compute_flip_normals_dim(Mesh* mesh, Bytes sides_are_exposed,
    Bytes verts_matter,
    Real simple_algorithm_threshold =
        0.95) {  // constant given by Aubry and Lohner
  constexpr auto side_dim = dim - 1;
  constexpr auto max_adj_sides =
      SimplexAvgDegree<side_dim, VERT, side_dim>::value * 4;
  auto v2s = mesh->ask_up(VERT, side_dim);
  auto sv2v = mesh->ask_verts_of(side_dim);
  auto verts_that_matter = collect_marked(verts_matter);
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * dim, 0.0);
  auto f = OMEGA_H_LAMBDA(LO vm) {
    auto v = verts_that_matter[vm];
    auto N_c = zero_vector<dim>();
    Few<Vector<dim>, max_adj_sides> N;
    auto n = 0;
    for (auto vs = v2s.a2ab[v]; vs < v2s.a2ab[v + 1]; ++vs) {
      auto s = v2s.ab2b[vs];
      if (!sides_are_exposed[s]) continue;
      OMEGA_H_CHECK(n < max_adj_sides);
      auto ssv2v = gather_verts<side_dim + 1>(sv2v, s);
      auto ssv2x = gather_vectors<side_dim + 1, dim>(coords, ssv2v);
      auto svec = get_side_vector(ssv2x);
      N_c += svec;  // svec is the area weighted face normal
      N[n++] = normalize(svec);
    }
    // as suggested by Aubry and Lohner, start with an initial guess done by
    // some kind of averaging (area weighted being an option)
    N_c = normalize(N_c);
    // if that is roughly the same as all the face normals, just accept it
    // and don't bother with the expensive algorithm.
    // this is especially a good idea because most models are mostly flat
    // surfaces
    Int i;
    for (i = 0; i < n; ++i) {
      if (N_c * N[i] < simple_algorithm_threshold) break;
    }
    if (i < n) {
      // nope, we actually have some nontrivial normals here.
      // run the super expensive algorithm.
      N_c = get_most_normal_normal(N, n);
    }
    set_vector(out, v, N_c);
  };
  parallel_for(verts_that_matter.size(), f, "compute_flip_normals");
  return mesh->sync_array(VERT, Reals(out), dim);
}

template <Int dim>
Bytes prevent_coarsen_flip2_dim(Mesh* mesh,
    Bytes sides_matter,  // exposed sides adjacent to candidate vertices
    Bytes verts_matter,  // vertices adjacent to sides that matter
    Reals vert_normals, LOs cands2edges, Bytes cand_codes,
    Real epsilon = OMEGA_H_EPSILON) {  // below this value a dot product is
                                       // considered negative
  OMEGA_H_CHECK(mesh->dim() == dim);
  constexpr auto side_dim = dim - 1;
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto sv2v = mesh->ask_verts_of(side_dim);
  auto v2s = mesh->ask_up(VERT, side_dim);
  auto coords = mesh->coords();
  auto ncands = cands2edges.size();
  auto out = Write<Byte>(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto e = cands2edges[cand];
    auto code = cand_codes[cand];
    for (Int eev_col = 0; eev_col < 2; ++eev_col) {
      if (!collapses(code, eev_col)) continue;
      auto v_col = ev2v[e * 2 + eev_col];
      if (!verts_matter[v_col]) continue;
      auto eev_onto = 1 - eev_col;
      auto v_onto = ev2v[e * 2 + eev_onto];
      for (auto vs = v2s.a2ab[v_col]; vs < v2s.a2ab[v_col + 1]; ++vs) {
        auto s = v2s.ab2b[vs];
        if (!sides_matter[s]) continue;
        auto vs_code = v2s.codes[vs];
        auto ssv_col = code_which_down(vs_code);
        auto ssv2v = gather_verts<side_dim + 1>(sv2v, s);
        Int ssv;
        for (ssv = 0; ssv < side_dim + 1; ++ssv) {
          if (ssv2v[ssv] == v_onto) break;  // ignore sides that will disappear
        }
        if (ssv != side_dim + 1)
          continue;               // ignore sides that will disappear (part 2)
        ssv2v[ssv_col] = v_onto;  // simulate the edge collapse on this side
        auto ssv2x = gather_vectors<side_dim + 1, dim>(coords, ssv2v);
        auto sn = get_side_vector(ssv2x);  // get its new normal
        // compare the new side normal with the old normals of its vertices
        auto ssv2n = gather_vectors<side_dim + 1, dim>(vert_normals, ssv2v);
        for (ssv = 0; ssv < side_dim + 1; ++ssv) {
          if (ssv2n[ssv] * sn < epsilon) break;
        }
        if (ssv < side_dim + 1) {  // got a violator here
          code = dont_collapse(code, eev_col);
        }
      }
    }
    out[cand] = code;
  };
  parallel_for(ncands, f, "prevent_coarsen_flip");
  return mesh->sync_subset_array(
      EDGE, Bytes(out), cands2edges, I8(DONT_COLLAPSE), 1);
}

template <Int dim>
Bytes prevent_coarsen_flip_dim(
    Mesh* mesh, LOs cands2edges, Bytes cand_codes) {
  auto edges_are_cands = mark_image(cands2edges, mesh->nedges());
  auto verts_are_cands = mark_down(mesh, EDGE, VERT, edges_are_cands);
  auto side_dim = dim - 1;
  auto sides_are_adj = mark_up(mesh, VERT, side_dim, verts_are_cands);
  auto sides_are_exposed = mark_exposed_sides(mesh);
  auto sides_matter = land_each(sides_are_adj, sides_are_exposed);
  auto verts_matter = mark_down(mesh, side_dim, VERT, sides_matter);
  auto vert_normals =
      compute_flip_normals_dim<dim>(mesh, sides_are_exposed, verts_matter);
  auto new_codes = prevent_coarsen_flip2_dim<dim>(
      mesh, sides_matter, verts_matter, vert_normals, cands2edges, cand_codes);
  return new_codes;
}

Bytes prevent_coarsen_flip(Mesh* mesh, LOs cands2edges, Bytes cand_codes) {
  if (mesh->dim() == 3) {
    return prevent_coarsen_flip_dim<3>(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 2) {
    return prevent_coarsen_flip_dim<2>(mesh, cands2edges, cand_codes);
  }
  return cand_codes;
}

}  // namespace Omega_h
