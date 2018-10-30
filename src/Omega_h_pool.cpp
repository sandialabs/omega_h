#include <Omega_h_pool.hpp>
#include <Omega_h_profile.hpp>
#include <Omega_h_fail.hpp>
#include <algorithm>

//DEBUG
#include <iostream>

namespace Omega_h {

static void call_underlying_frees(Pool& pool, BlockList& list) {
  for (auto& block : list) {
    pool.underlying_free(block.data);
  }
  list.clear();
}

Pool::Pool(std::size_t page_size_in,
    MallocFunc malloc_in,
    FreeFunc free_in)
  :page_size(page_size_in)
  ,underlying_malloc(malloc_in)
  ,underlying_free(free_in)
{
}

Pool::~Pool() {
  call_underlying_frees(*this, used_blocks);
  call_underlying_frees(*this, free_blocks);
}

static BlockList::iterator find_best_fit(Pool& pool, std::size_t size) {
  auto const end = pool.free_blocks.end();
  auto best = end;
  for (auto it = pool.free_blocks.begin(); it != end; ++it) {
    if (it->size < size) continue;
    if ((it->size > pool.page_size) && ((it->size / 2) > size)) continue;
    if (best == end || it->size < best->size) {
      best = it;
    }
  }
  return best;
}

static std::size_t underlying_total_size(Pool& pool) {
  std::size_t total_size = 0;
  for (auto& block : pool.free_blocks) {
    total_size += block.size;
  }
  for (auto& block : pool.used_blocks) {
    total_size += block.size;
  }
  return total_size;
}

void* allocate(Pool& pool, std::size_t size) {
  ScopedTimer timer("pool allocate");
  auto const best_fit = find_best_fit(pool, size);
  if (best_fit != pool.free_blocks.end()) {
    auto const data = best_fit->data;
    std::cerr << "pool resurrected " << data << ", " << best_fit->size << " for " << size << '\n';
    pool.used_blocks.push_back(*best_fit);
    pool.free_blocks.erase(best_fit);
    return data;
  }
  auto pages = size / pool.page_size;
  if (size % pool.page_size) pages += 1;
  auto const size_to_alloc = pages * pool.page_size;
  auto data = pool.underlying_malloc(size_to_alloc);
  if (data == nullptr) {
    call_underlying_frees(pool, pool.free_blocks);
    data = pool.underlying_malloc(size_to_alloc);
  }
  if (data == nullptr) {
    Omega_h_fail("Pool failed to allocate %zu bytes, %zu bytes already allocated\n",
        size_to_alloc, underlying_total_size(pool));
  }
  pool.used_blocks.push_back({size_to_alloc, data});
  std::cerr << "pool allocated " << data << ", " << size_to_alloc << " for " << size << '\n';
  return data;
}

void deallocate(Pool& pool, void* data) {
  ScopedTimer timer("pool deallocate");
  auto const end = pool.used_blocks.end();
  auto const it = std::find_if(pool.used_blocks.begin(), end, [=](Block const& b) -> bool { return b.data == data; });
  std::cerr << "pool recycling " << data << ", " << it->size << '\n';
  if (it == end) {
    Omega_h_fail("Tried to deallocate %p from pool, but pool didn't allocate it\n", data);
  }
  pool.free_blocks.push_back(*it);
  pool.used_blocks.erase(it);
}

}
