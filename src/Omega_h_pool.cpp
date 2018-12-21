#include <Omega_h_fail.hpp>
#include <Omega_h_pool.hpp>
#include <Omega_h_profile.hpp>
#include <algorithm>

namespace Omega_h {

static void call_underlying_frees(Pool& pool, BlockList list[]) {
  for (std::size_t i = 0; i < 64; ++i) {
    for (auto block : list[i]) {
      pool.underlying_free(block, (std::size_t(1) << i));
    }
    list[i].clear();
  }
}

Pool::Pool(MallocFunc malloc_in, FreeFunc free_in)
    : underlying_malloc(malloc_in), underlying_free(free_in) {}

Pool::~Pool() {
  call_underlying_frees(*this, used_blocks);
  call_underlying_frees(*this, free_blocks);
}

static std::size_t underlying_total_size(Pool& pool) {
  std::size_t total_size = 0;
  for (std::size_t i = 0; i < 64; ++i) {
    total_size += (std::size_t(1) << i) * pool.used_blocks[i].size();
    total_size += (std::size_t(1) << i) * pool.free_blocks[i].size();
  }
  return total_size;
}

void* allocate(Pool& pool, std::size_t size) {
  ScopedTimer timer("pool allocate");
  std::size_t shift;
  for (shift = 0; ((std::size_t(1) << shift) < size); ++shift)
    ;
  if (!pool.free_blocks[shift].empty()) {
    auto const data = pool.free_blocks[shift].back();
    pool.used_blocks[shift].push_back(data);
    pool.free_blocks[shift].pop_back();
    return data;
  }
  auto const size_to_alloc = (std::size_t(1) << shift);
  auto data = pool.underlying_malloc(size_to_alloc);
  if (data == nullptr) {
    call_underlying_frees(pool, pool.free_blocks);
    data = pool.underlying_malloc(size_to_alloc);
  }
  if (data == nullptr) {
    Omega_h_fail(
        "Pool failed to allocate %zu bytes, %zu bytes already allocated\n",
        size_to_alloc, underlying_total_size(pool));
  }
  pool.used_blocks[shift].push_back(data);
  return data;
}

void deallocate(Pool& pool, void* data, std::size_t size) {
  ScopedTimer timer("pool deallocate");
  std::size_t shift;
  for (shift = 0; ((std::size_t(1) << shift) < size); ++shift)
    ;
  auto const end = pool.used_blocks[shift].end();
  auto const it = std::find_if(pool.used_blocks[shift].begin(), end,
      [=](VoidPtr b) -> bool { return b == data; });
  if (it == end) {
    Omega_h_fail(
        "Tried to deallocate %p from pool, but pool didn't allocate it\n",
        data);
  }
  pool.free_blocks[shift].push_back(*it);
  pool.used_blocks[shift].erase(it);
}
}  // namespace Omega_h
