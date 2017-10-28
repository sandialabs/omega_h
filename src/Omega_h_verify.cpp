#include <Omega_h_verify.hpp>

#include <Omega_h_align.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_simplex.hpp>
#include <Omega_h_loop.hpp>

namespace Omega_h {

OMEGA_H_INLINE Int down_template2(
    Int high_dim, Int mid_dim, Int low_dim, Int which_mid, Int which_low) {
  switch (high_dim) {
    case 1:
      switch (mid_dim) {
        case 0:
          return which_mid;
      }
    case 2:
      switch (mid_dim) {
        case 0:
          return which_mid;
        case 1:
          switch (low_dim) {
            case 0:
              switch (which_mid) {
                case 0:
                  switch (which_low) {
                    case 0:
                      return 0;
                    case 1:
                      return 1;
                  }
                case 1:
                  switch (which_low) {
                    case 0:
                      return 1;
                    case 1:
                      return 2;
                  }
                case 2:
                  switch (which_low) {
                    case 0:
                      return 2;
                    case 1:
                      return 0;
                  }
              }
            case 1:
              return which_mid;
          }
      }
    case 3:
      switch (mid_dim) {
        case 0:
          return which_mid;
        case 1:
          switch (low_dim) {
            case 0:
              switch (which_mid) {
                case 0:
                  switch (which_low) {
                    case 0:
                      return 0;
                    case 1:
                      return 1;
                  }
                case 1:
                  switch (which_low) {
                    case 0:
                      return 1;
                    case 1:
                      return 2;
                  }
                case 2:
                  switch (which_low) {
                    case 0:
                      return 2;
                    case 1:
                      return 0;
                  }
                case 3:
                  switch (which_low) {
                    case 0:
                      return 0;
                    case 1:
                      return 3;
                  }
                case 4:
                  switch (which_low) {
                    case 0:
                      return 1;
                    case 1:
                      return 3;
                  }
                case 5:
                  switch (which_low) {
                    case 0:
                      return 2;
                    case 1:
                      return 3;
                  }
              }
            case 1:
              return which_mid;
          }
        case 2:
          switch (low_dim) {
            case 0:
              switch (which_mid) {
                case 0:
                  switch (which_low) {
                    case 0:
                      return 0;
                    case 1:
                      return 2;
                    case 2:
                      return 1;
                  }
                case 1:
                  switch (which_low) {
                    case 0:
                      return 0;
                    case 1:
                      return 1;
                    case 2:
                      return 3;
                  }
                case 2:
                  switch (which_low) {
                    case 0:
                      return 1;
                    case 1:
                      return 2;
                    case 2:
                      return 3;
                  }
                case 3:
                  switch (which_low) {
                    case 0:
                      return 2;
                    case 1:
                      return 0;
                    case 2:
                      return 3;
                  }
              }
            case 1:
              switch (which_mid) {
                case 0:
                  switch (which_low) {
                    case 0:
                      return 2;
                    case 1:
                      return 1;
                    case 2:
                      return 0;
                  }
                case 1:
                  switch (which_low) {
                    case 0:
                      return 0;
                    case 1:
                      return 4;
                    case 2:
                      return 3;
                  }
                case 2:
                  switch (which_low) {
                    case 0:
                      return 1;
                    case 1:
                      return 5;
                    case 2:
                      return 4;
                  }
                case 3:
                  switch (which_low) {
                    case 0:
                      return 2;
                    case 1:
                      return 3;
                    case 2:
                      return 5;
                  }
              }
          }
      }
  }
  return -1;
}

/* verify that the downward vertex arrays are
   consistent with one another.
   e.g. the first face of tetrahedron {4,8,7,6}
   must be {4,8,7}. */
void verify_down_verts(Mesh* mesh) {
  for (Int ld = 0; ld < 1; ++ld) {
    for (Int md = ld + 1; md < mesh->dim(); ++md) {
      for (Int hd = md + 1; hd <= mesh->dim(); ++hd) {
        auto h2l =  mesh->ask_down(hd, ld);
        auto h2m =  mesh->ask_down(hd, md);
        auto m2l =  mesh->ask_down(md, ld);
        auto nhhm = simplex_degree(hd, md);
        auto nhhl = simplex_degree(hd, ld);
        auto nmml = simplex_degree(md, ld);
        auto f = OMEGA_H_LAMBDA(LO h) {
          for (Int hhm = 0; hhm < nhhm; ++hhm) {
            auto m = h2m.ab2b[h * nhhm + hhm];
            auto code = h2m.codes[h * nhhm + hhm];
            for (Int mml = 0; mml < nmml; ++mml) {
              auto l = m2l.ab2b[m * nmml + mml];
              auto hml = align_index(nmml, ld, mml, code);
              auto hhl = down_template2(hd, md, ld, hhm, hml);
              auto l2 = h2l.ab2b[h * nhhl + hhl];
              OMEGA_H_CHECK(l == l2);
            }
          }
        };
        parallel_for(mesh->nents(hd), f);
      }
    }
  }
}

}  // end namespace Omega_h
