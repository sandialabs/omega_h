#include <Omega_h_verify.hpp>

#include <Omega_h_align.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_simplex.hpp>
#include <Omega_h_loop.hpp>
#include <Omega_h_adj.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_array_ops.hpp>

#include <iostream>

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
bool verify_down_verts(Mesh* mesh) {
  for (Int ld = 0; ld <= 1; ++ld) {
    for (Int md = ld + 1; md < mesh->dim(); ++md) {
      for (Int hd = md + 1; hd <= mesh->dim(); ++hd) {
      //std::cerr << "verify_down_verts ld " << ld << " md " << md << " hd " << hd << '\n';
        auto h2l =  mesh->ask_down(hd, ld);
        auto h2m =  mesh->ask_down(hd, md);
        auto m2l =  mesh->ask_down(md, ld);
        auto nhhm = simplex_degree(hd, md);
        auto nhhl = simplex_degree(hd, ld);
        auto nmml = simplex_degree(md, ld);
      //auto f = OMEGA_H_LAMBDA(LO h) {
        for (LO h = 0; h < mesh->nents(hd); ++h) {
          for (Int hhm = 0; hhm < nhhm; ++hhm) {
            auto m = h2m.ab2b[h * nhhm + hhm];
          //if (ld == 1 && h == 143137) {
          //  std::cerr << "h2m.ab2b[" << h << " * " << nhhm << " + " << hhm << "] = " << m << '\n';
          //}
            auto code = h2m.codes[h * nhhm + hhm];
            for (Int mml = 0; mml < nmml; ++mml) {
              auto l = m2l.ab2b[m * nmml + mml];
            //if (ld == 1 && h == 143137) {
            //  std::cerr << "m2l.ab2b[" << m << " * " << nmml << " + " << mml << "] = " << l << '\n';
            //}
              auto hml = align_index(nmml, ld, mml, code);
              auto hhl = down_template2(hd, md, ld, hhm, hml);
              auto l2 = h2l.ab2b[h * nhhl + hhl];
            //if (ld == 1 && h == 143137) {
            //  std::cerr << "h2l.ab2b[" << h << " * " << nhhl << " + " << hhl << "] = " << l2 << '\n';
            //}
              if (l != l2) {
                std::cerr << "hd " << hd << " md " << md << " ld " << ld << '\n';
                std::cerr << "h " << h << " m " << m << " l " << l << '\n';
                auto edges2verts = mesh->ask_down(1, 0).ab2b;
                std::cerr << "local edge " << hhl << " of element " << h << " is " << l2
                  << " which connects " << edges2verts[l2 * 2 + 0]
                  << " and " << edges2verts[l2 * 2 + 1] << '\n';
                std::cerr << "local edge " << mml << " of face " << m << " is " << l
                  << " which connects " << edges2verts[l * 2 + 0]
                  << " and " << edges2verts[l * 2 + 1] << '\n';
                return false;
              }
              OMEGA_H_CHECK(l == l2);
            }
          }
        }
      //};
      //parallel_for(mesh->nents(hd), f);
      }
    }
  }
  std::cout << "mesh verified!\n";
  return true;
}

void verify_no_duplicates(Mesh* mesh) {
  for (Int ent_dim = 1; ent_dim <= mesh->dim(); ++ent_dim) {
    std::cerr << "checking for duplicates in dim " << ent_dim << " entities...\n";
    auto ev2v = mesh->ask_verts_of(ent_dim);
    auto v2e = mesh->ask_up(VERT, ent_dim);
    LOs a2b;
    Bytes codes;
    find_matches(ent_dim, ev2v, ev2v, v2e, &a2b, &codes);
    std::cerr << "no duplicates in dim " << ent_dim << " entities...\n";
  }
}

void verify_class(Mesh* mesh) {
  for (Int ent_dim = mesh->dim() - 1; ent_dim >= VERT; --ent_dim) {
    std::cerr << "verifying classification, dim " << ent_dim << '\n';
    auto class_dims = mesh->get_array<Byte>(ent_dim, "class_dim");
    auto class_ids = mesh->get_array<ClassId>(ent_dim, "class_id");
    auto class_dims_w = deep_copy(class_dims);
    auto class_ids_w = deep_copy(class_ids);
    project_classification(mesh, ent_dim, class_dims_w, class_ids_w);
    if (ent_dim == 1) {
      auto ev2v = mesh->ask_verts_of(ent_dim);
      auto vert_class_dims = mesh->get_array<Byte>(VERT, "class_dim");
      auto vert_class_ids = mesh->get_array<ClassId>(VERT, "class_id");
      auto tri_class_dims = mesh->get_array<Byte>(TRI, "class_dim");
      auto tri_class_ids = mesh->get_array<ClassId>(TRI, "class_id");
      auto tet_class_dims = mesh->get_array<Byte>(TET, "class_dim");
      auto tet_class_ids = mesh->get_array<ClassId>(TET, "class_id");
      auto e2t = mesh->ask_up(EDGE, TRI);
      auto t2k = mesh->ask_up(TRI, TET);
      for (LO i = 0; i < mesh->nents(ent_dim); ++i) {
        if ((ev2v[i * 2 + 0] == 188224 && ev2v[i * 2 + 1] == 293364) ||
            (ev2v[i * 2 + 0] == 293364 && ev2v[i * 2 + 1] == 188224) ||
            (ev2v[i * 2 + 0] == 293364 && ev2v[i * 2 + 1] == 19882) ||
            (ev2v[i * 2 + 0] == 19882 && ev2v[i * 2 + 1] == 293364)) {
          std::cerr << "edge " << i << " has vertices "
            << ev2v[i * 2 + 0] << " and " << ev2v[i * 2 + 1] << '\n';
          std::cerr << "edge " << i << " is classified on model " << class_ids[i]
            << " of dimension " << Int(class_dims[i]) << '\n';
          for (Int j = 0; j < 2; ++j) {
            auto v = ev2v[i * 2 + j];
            std::cerr << "vertex " << v << " is classified on model " << vert_class_ids[v]
              << " of dimension " << Int(vert_class_dims[v]) << '\n';
          }
          for (auto et = e2t.a2ab[i]; et < e2t.a2ab[i + 1]; ++et) {
            auto t = e2t.ab2b[et];
            std::cerr << "edge " << i << " is adjacent to triangle " << t
              << " which is classified on model " << tri_class_ids[t]
              << " of dimension " << Int(tri_class_dims[t]) << '\n';
            for (auto tk = t2k.a2ab[t]; tk < t2k.a2ab[t + 1]; ++tk) {
              auto k = t2k.ab2b[tk];
              std::cerr << "triangle " << t << " is adjacent to tetrahedron "
                << k << " which is classified on model " << tet_class_ids[k]
                << " of dimension " << Int(tet_class_dims[k]) << '\n';
            }
          }
        }
      }
    }
    if (get_min(each_eq(class_dims, Read<Byte>(class_dims_w))) == Byte(0)) {
      Omega_h_fail("class_dim for dimension %d is inconsistent\n", ent_dim);
    }
    if (get_min(each_eq(class_ids, Read<ClassId>(class_ids_w))) == Byte(0)) {
      Omega_h_fail("class_id for dimension %d is inconsistent\n", ent_dim);
    }
  }
}

}  // end namespace Omega_h
