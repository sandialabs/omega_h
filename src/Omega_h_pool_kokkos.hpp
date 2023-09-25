//
// Created by Matthew McCall on 6/12/23.
//
// Derived from
// https://github.com/matthew-mccall/kokkos-memory-pool/blob/8e0a45a5b5d6823976d867e20d742314b6a1d268/src/MemoryPool/MemoryPool.hpp

#ifndef OMEGA_H_KOKKOS_POOL_HPP
#define OMEGA_H_KOKKOS_POOL_HPP

#include <cstddef>
#include <list>
#include <map>
#include <set>
#include <utility>

#include "Kokkos_Core.hpp"
#include "Omega_h_kokkos.hpp"

namespace Omega_h {

constexpr size_t DEFAULT_BYTES_PER_CHUNK = 1000;

// Total initial size in bytes will subsequently be DEFAULT_NUMBER_OF_CHUNKS * DEFAULT_BYTES_PER_CHUNK
constexpr size_t DEFAULT_NUMBER_OF_CHUNKS = 700'000;

// If the largest StaticKokkosPool is 2Kb, then the next StaticKokkosPool will be at least 4kB resulting in a total pool size of 2Kb + 4Kb = 6Kb
constexpr size_t DEFAULT_GROWTH_CONSTANT = 2;

using IndexPair = std::pair<size_t, size_t>;

class CompareFreeIndices {
 public:
  using is_transparent =
      void;  // https://www.fluentcpp.com/2017/06/09/search-set-another-type-key/

  auto operator()(IndexPair lhs, IndexPair rhs) const -> bool;
  auto operator()(IndexPair lhs, size_t rhs) const -> bool;
  auto operator()(size_t lhs, IndexPair rhs) const -> bool;
};

using MultiSetBySizeT = std::multiset<IndexPair, CompareFreeIndices>;
using SetByIndexT = std::set<IndexPair>;

/**
 * @brief A memory pool for allocating and deallocating chunks of memory.
 *
 * Does not support resizing.
 */
class StaticKokkosPool {
 public:
  StaticKokkosPool(size_t numChunks, size_t bytesPerChunk);

  auto allocate(size_t n) -> void*;
  void deallocate(void* data);

  auto getNumAllocations() const -> unsigned;
  auto getNumFreeChunks() const -> unsigned;
  auto getNumAllocatedChunks() const -> unsigned;
  auto getNumChunks() const -> unsigned;
  auto getNumFreeFragments() const -> unsigned;

  static auto getRequiredChunks(size_t n, size_t bytesPerChunk) -> size_t;

  ~StaticKokkosPool();

 private:
  auto insertIntoSets(IndexPair indices)
      -> std::pair<MultiSetBySizeT::iterator, SetByIndexT::iterator>;
  void removeFromSets(IndexPair indices);

  const size_t numberOfChunks;
  const size_t chunkSize;

  void* pool;
  MultiSetBySizeT freeSetBySize;  // For finding free chunks logarithmically
  SetByIndexT freeSetByIndex;     // For merging adjacent free chunks
  std::map<void*, IndexPair> allocations;
};

/**
 * @brief A memory pool for allocating and deallocating chunks of memory.
 *
 * It is a collection of StaticKokkosPool objects. This allows for resizing.
 */
class KokkosPool {
 public:
  explicit KokkosPool(size_t bytesPerChunks);

  auto allocate(size_t n) -> void*;
  void deallocate(void* data);

  template <typename DataType>
  auto allocateView(size_t n) -> View<DataType*> {
    return View<DataType*>(
        reinterpret_cast<DataType*>(allocate(n * sizeof(DataType))), n);
  }

  template <typename DataType>
  void deallocateView(View<DataType*> view) {
    deallocate(reinterpret_cast<void*>(view.data()));
  }

  auto getNumAllocations() const -> unsigned;
  auto getNumFreeChunks() const -> unsigned;
  auto getNumAllocatedChunks() const -> unsigned;
  auto getNumChunks() const -> unsigned;
  auto getNumFreeFragments() const -> unsigned;
  auto getChunkSize() const -> size_t;

  static auto getGlobalPool() -> KokkosPool&;
  static auto destroyGlobalPool() -> void;

 private:
  using PoolListT = std::list<StaticKokkosPool>;

  size_t chunkSize;
  PoolListT pools;
  std::map<void*, PoolListT::iterator> allocations;
};

}  // namespace Omega_h

#endif  // OMEGA_H_KOKKOS_POOL_HPP
