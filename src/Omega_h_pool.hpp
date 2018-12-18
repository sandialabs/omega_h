#ifndef OMEGA_H_POOL_HPP
#define OMEGA_H_POOL_HPP

#include <functional>
#include <vector>

namespace Omega_h {

using VoidPtr = void*;
using BlockList = std::vector<VoidPtr>;
using MallocFunc = std::function<VoidPtr(std::size_t)>;
using FreeFunc = std::function<void(VoidPtr, std::size_t)>;

struct Pool {
  Pool(MallocFunc, FreeFunc);
  ~Pool();
  Pool(Pool const&) = delete;
  Pool(Pool&&) = delete;
  Pool& operator=(Pool const&) = delete;
  Pool& operator=(Pool&&) = delete;
  BlockList used_blocks[64];
  BlockList free_blocks[64];
  MallocFunc underlying_malloc;
  FreeFunc underlying_free;
};

void* allocate(Pool&, std::size_t);
void deallocate(Pool&, void*, std::size_t);
}  // namespace Omega_h

#endif
