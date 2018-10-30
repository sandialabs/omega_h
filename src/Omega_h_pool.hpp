#ifndef OMEGA_H_POOL_HPP
#define OMEGA_H_POOL_HPP

#include <vector>
#include <functional>

namespace Omega_h {

struct Block {
  std::size_t size;
  void* data;
};

using BlockList = std::vector<Block>;
using VoidPtr = void*;
using MallocFunc = std::function<VoidPtr(std::size_t)>;
using FreeFunc = std::function<void(VoidPtr)>;

struct Pool {
  Pool(std::size_t, MallocFunc, FreeFunc);
  ~Pool();
  Pool(Pool const&) = delete;
  Pool(Pool&&) = delete;
  Pool& operator=(Pool const&) = delete;
  Pool& operator=(Pool&&) = delete;
  std::size_t page_size;
  BlockList used_blocks;
  BlockList free_blocks;
  MallocFunc underlying_malloc;
  FreeFunc underlying_free;
};

void* allocate(Pool&, std::size_t);
void deallocate(Pool&, void*);

}

#endif
