#include <iostream>

#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_timer.hpp>


#include <Omega_h_for.hpp>

void call_print(Omega_h::LOs a) {
  printf("\n");
  auto a_w = Omega_h::Write<Omega_h::LO> (a.size());
  auto r2w = OMEGA_H_LAMBDA(Omega_h::LO i) {
    a_w[i] = a[i];
  };
  Omega_h::parallel_for(a.size(), r2w);
  auto a_host = Omega_h::HostWrite<Omega_h::LO>(a_w);
  for (int i=0; i<a_host.size(); ++i) {
    printf(" %d,", a_host[i]);
  };
  printf("\n");
  printf("\n");
  return;
}

void print_owners(Omega_h::Remotes owners, int rank) {
  printf("\n");
  auto ranks = owners.ranks;
  auto idxs = owners.idxs;
  auto ranks_w = Omega_h::Write<Omega_h::LO> (ranks.size());
  auto idxs_w = Omega_h::Write<Omega_h::LO> (idxs.size());
  auto r2w = OMEGA_H_LAMBDA(Omega_h::LO i) {
    ranks_w[i] = ranks[i];
    idxs_w[i] = idxs[i];
  };  
  Omega_h::parallel_for(idxs.size(), r2w);
  auto ranks_host = Omega_h::HostWrite<Omega_h::LO>(ranks_w);
  auto idxs_host = Omega_h::HostWrite<Omega_h::LO>(idxs_w);
  printf("On rank %d\n", rank);
  for (int i=0; i<idxs_host.size(); ++i) {
    printf("owner of %d, is on rank %d, with LId %d\n", i, ranks_host[i], idxs_host[i]);
  };  
  printf("\n");
  printf("\n");
  return;
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  if (argc != 5) {
    if (!world->rank()) {
      std::cout << "usage: " << argv[0] << " in.osh <nparts> out.osh model_file\n";
    }
    return -1;
  }
  auto nparts_total = world->size();
  auto path_in = argv[1];
  auto nparts_out = atoi(argv[2]);
  auto path_out = argv[3];
  auto path_model = argv[4];
  auto t0 = Omega_h::now();
  if (nparts_out < 1) {
    if (!world->rank()) {
      std::cout << "error: invalid output part count " << nparts_out << '\n';
    }
    return -1;
  }
  if (nparts_out > nparts_total) {
    if (!world->rank()) {
      std::cout << "error: output part count " << nparts_out
                << " greater than MPI job size " << nparts_total << '\n';
    }
    return -1;
  }
  auto nparts_in = Omega_h::binary::read_nparts(path_in, world);
  if (nparts_in > nparts_total) {
    if (!world->rank()) {
      std::cout << "error: input part count " << nparts_in
                << " greater than MPI job size " << nparts_total << '\n';
    }
    return -1;
  }
  auto is_in = (world->rank() < nparts_in);
  auto comm_in = world->split(int(is_in), 0);
  auto is_out = (world->rank() < nparts_out);
  auto comm_out = world->split(int(is_out), 0);
  auto mesh = Omega_h::Mesh(&lib);
  auto version = Omega_h::binary::read_version(path_in, world);
  if (is_in) {
    Omega_h::binary::read_in_comm(path_in, comm_in, &mesh, version);
    Omega_h::binary::read_model(path_model, &mesh);
    if (nparts_out < nparts_in) {
      Omega_h_fail(
          "partitioning to a smaller part count not yet implemented\n");
    }
  }
/*
  if (!world->rank()) {
    for (int d=0; d<mesh.dim(); ++d) {
      printf("dim=%d\n", d);
      printf("owners\n");
      print_owners(mesh.ask_owners(d), world->rank());
      printf("model_ents\n");
      call_print(mesh.ask_model_ents(d));
      printf("model_matches\n");
      call_print(mesh.ask_model_matches(d));
    }
  }
*/
  //printf("ok1\n");
  if (is_in || is_out) mesh.set_comm(comm_out); //dist changes here but owners same as before
/*
  if (!world->rank()) {
    for (int d=0; d<mesh.dim(); ++d) {
      printf("dim=%d\n", d);
      printf("owners with new dist\n");
      print_owners(mesh.ask_owners(d), world->rank());
    }
  }
*/
  //printf("ok2\n");
  //int waiting=1;
  //while(waiting);
  if (is_out) {
    if (nparts_out != nparts_in) mesh.balance();
    Omega_h::binary::write(path_out, &mesh);
  }
  world->barrier();
  auto t1 = Omega_h::now();
  auto imb = mesh.imbalance();
  if (!world->rank()) {
    std::cout << "repartitioning took " << (t1 - t0) << " seconds\n";
    std::cout << "imbalance is " << imb << "\n";
  }
}
