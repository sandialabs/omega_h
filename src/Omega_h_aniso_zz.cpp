#include "Omega_h_size_field.hpp"

#include "access.hpp"
#include "algebra.hpp"
#include "derive.hpp"
#include "simplices.hpp"
#include "Omega_h_few.hpp"
#include "size.hpp"

/* Micheletti, S., and S. Perotto.
 * "Anisotropic adaptation via a Zienkiewicz–Zhu error estimator
 *  for 2D elliptic problems."
 * Numerical Mathematics and Advanced Applications 2009.
 * Springer Berlin Heidelberg, 2010. 645-653.
 *
 * Farrell, P. E., S. Micheletti, and S. Perotto.
 * "An anisotropic Zienkiewicz–Zhu‐type error estimator for 3D applications."
 * International journal for numerical methods in engineering
 * 85.6 (2011): 671-692.
 */

template <Int dim>
Reals get_aniso_zz_metric_dim(Mesh* mesh, Reals elem_gradients,
    Real error_bound, Real max_size) {
  constexpr auto nverts_per_elem = dim + 1;
  auto elem_verts2verts = mesh->ask_elem_verts();
  auto verts2elems = mesh->ask_up(VERT, dim);
  constexpr auto max_elems_per_patch =
    AvgDegree<dim, 0, dim>::value * nverts_per_elem;
  auto elems2volume = measure_elements_real(mesh);
  auto nglobal_elems = get_sum(mesh->comm(), mesh->owned(dim));
  auto iso_inv = get_iso_simplex_inv<dim>();
  auto coords = mesh->coords();
  auto f = LAMBDA(LO elem) {
    auto k_verts = gather_verts<dim + 1>(elem_verts2vert, elem);
    auto p = gather_vectors<dim + 1, dim>(coords, k_verts);
    auto m_k = iso_inv * get_centered_basis(p);
    Few<LO, max_elems_per_patch> patch_elems;
    Int npatch_elems = 0;
    for (auto ev = elem * nverts_per_elem;
         ev < ((elem + 1) * nverts_per_elem); ++ev) {
      auto vert = elem_verts2vert[ev];
      for (auto ve = verts2elems.a2ab[vert];
           ve < verts2elems.a2ab[vert + 1]; ++ve) {
        auto patch_elem = verts2elems.ab2b[ve];
        add_unique(patch_elems, npatch_elems, patch_elem);
      }
    }
    Real patch_volume = 0;
    auto op_sum = zero_matrix<dim, dim>();
    for (Int i = 0; i < npatch_elems; ++i) {
      auto patch_elem = patch_elems[i];
      auto gradient = get_vector<dim>(elem_gradients, patch_elem);
      auto volume = elems2volume[patch_elem];
      auto op = outer_product(gradient, gradient);
      patch_volume += volume;
      op_sum = op_sum + (op * volume);
    }
    auto patch_op = op_sum / patch_volume;
    auto decomp = decompose_eigen(patch_op);
    auto g = decomp.l;
    auto a = square(error_bound) / (dim * nglobal_elems);
  };
}
