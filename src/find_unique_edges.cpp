#include "Omega_h_swap3d_tables.hpp"

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

using namespace Omega_h::swap3d;

static int unique_push_back(
    std::vector<std::pair<int, int>>& v, std::pair<int, int> p) {
  auto it = std::find(v.begin(), v.end(), p);
  if (it == v.end()) {
    v.push_back(p);
    return int(v.size() - 1);
  } else {
    return int(it - v.begin());
  }
}

int main() {
  std::vector<std::vector<std::pair<int, int>>> unique_edges(4);
  std::vector<std::vector<std::vector<int>>> edges2unique_edges(4);
  for (int loop_size = 4; loop_size <= MAX_EDGE_SWAP; ++loop_size) {
    unique_edges.push_back(std::vector<std::pair<int, int>>());
    edges2unique_edges.push_back(std::vector<std::vector<int>>());
    auto nedges = swap_nint_edges[loop_size];
    auto nmeshes = swap_mesh_counts[loop_size];
    for (int mesh = 0; mesh < nmeshes; ++mesh) {
      edges2unique_edges.back().push_back(std::vector<int>());
      for (int edge = 0; edge < nedges; ++edge) {
        auto v0 = swap_int_edges[loop_size][mesh][edge * 2 + 0];
        auto v1 = swap_int_edges[loop_size][mesh][edge * 2 + 1];
        auto ue =
            unique_push_back(unique_edges.back(), std::pair<int, int>(v0, v1));
        edges2unique_edges.back().back().push_back(ue);
      }
    }
  }
  std::cout << "typedef int IntPair[2];\n\n";
  for (std::size_t loop_size = 4; loop_size <= MAX_EDGE_SWAP; ++loop_size) {
    std::cout << "CONSTANT static IntPair const unique_edges_" << loop_size
              << "[" << unique_edges[loop_size].size() << "] = {\n";
    for (std::size_t edge = 0; edge < unique_edges[loop_size].size(); ++edge) {
      if (edge != 0) std::cout << ", ";
      std::cout << "{" << unique_edges[loop_size][edge].first << ", ";
      std::cout << unique_edges[loop_size][edge].second << "}\n";
    }
    std::cout << "};\n\n";
  }
  std::cout << "CONSTANT static IntPair const* const unique_edges["
            << "MAX_EDGE_SWAP + 1] = {\n";
  for (int loop_size = 0; loop_size < 4; ++loop_size) {
    if (loop_size != 0) std::cout << ", ";
    std::cout << "0" << '\n';
  }
  for (int loop_size = 4; loop_size <= MAX_EDGE_SWAP; ++loop_size) {
    std::cout << ", ";
    std::cout << "unique_edges_" << loop_size << '\n';
  }
  std::cout << "};\n\n";
  for (std::size_t loop_size = 4; loop_size <= MAX_EDGE_SWAP; ++loop_size) {
    auto nedges = std::size_t(swap_nint_edges[loop_size]);
    auto nmeshes = std::size_t(swap_mesh_counts[loop_size]);
    for (std::size_t mesh = 0; mesh < nmeshes; ++mesh) {
      std::cout << "CONSTANT static Int const edges2unique_" << loop_size << "_"
                << mesh << "[" << nedges << "] = {\n";
      for (std::size_t edge = 0; edge < nedges; ++edge) {
        if (edge != 0) std::cout << ", ";
        std::cout << edges2unique_edges[loop_size][mesh][edge] << '\n';
      }
      std::cout << "};\n\n";
    }
  }
  for (int loop_size = 4; loop_size <= MAX_EDGE_SWAP; ++loop_size) {
    auto nmeshes = swap_mesh_counts[loop_size];
    std::cout << "CONSTANT static Int const* const edges2unique_" << loop_size
              << "[" << nmeshes << "] = {\n";
    for (int mesh = 0; mesh < nmeshes; ++mesh) {
      if (mesh != 0) std::cout << ", ";
      std::cout << "edges2unique_" << loop_size << "_" << mesh << "\n";
    }
    std::cout << "};\n\n";
  }
  std::cout << "CONSTANT static Int const* const* const edges2unique["
            << "MAX_EDGE_SWAP + 1] = {\n";
  for (int loop_size = 0; loop_size < 4; ++loop_size) {
    if (loop_size != 0) std::cout << ", ";
    std::cout << "0" << '\n';
  }
  for (int loop_size = 4; loop_size <= MAX_EDGE_SWAP; ++loop_size) {
    std::cout << ", ";
    std::cout << "edges2unique_" << loop_size << '\n';
  }
  std::cout << "};\n\n";
}
