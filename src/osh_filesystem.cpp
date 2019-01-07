#include <Omega_h_filesystem.hpp>
#include <iostream>
#include <cstring>

int main(int argc, char** argv) {
  if (argc == 2 && 0 == strcmp(argv[1], "pwd")) {
    std::cout << Omega_h::filesystem::current_path() << '\n';
  }
  else if (argc == 2 && 0 == strcmp(argv[1], "ls")) {
    for (Omega_h::filesystem::directory_iterator it(Omega_h::filesystem::current_path()), end;
        it != end; ++it) {
      std::cout << it->path() << '\n';
    }
  }
  else if (argc == 3 && 0 == strcmp(argv[1], "rm")) {
    Omega_h::filesystem::remove(argv[2]);
  }
  else {
    std::cout << "usage: osh_filesystem pwd\n";
    std::cout << "       osh_filesystem ls\n";
    std::cout << "       osh_filesystem rm <file>\n";
  }
}
