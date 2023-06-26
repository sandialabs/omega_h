//
// Created by Matthew McCall on 6/12/23.
//
// Derived from
// https://github.com/matthew-mccall/kokkos-memory-pool/blob/8e0a45a5b5d6823976d867e20d742314b6a1d268/src/MemoryPool/MemoryPool.cpp

#include "Omega_h_pool_kokkos.hpp"

#include <cassert>
#include <cmath>
#include <iterator>
#include <optional>
#include <vector>

namespace {
std::optional<Omega_h::KokkosPool> s_pool;
}

namespace Omega_h {

auto CompareFreeIndices::operator()(IndexPair lhs, IndexPair rhs) const
    -> bool {
  auto [lhsStart, lhsEnd] = lhs;
  auto [rhsStart, rhsEnd] = rhs;
  size_t lhsSize = lhsEnd - lhsStart;
  size_t rhsSize = rhsEnd - rhsStart;

  if (lhsSize == rhsSize) {
    return lhs < rhs;
  }

  return lhsSize < rhsSize;
}

auto CompareFreeIndices::operator()(IndexPair lhs, size_t rhs) const -> bool {
  auto [lhsStart, lhsEnd] = lhs;
  return (lhsEnd - lhsStart) < rhs;
}

auto CompareFreeIndices::operator()(size_t lhs, IndexPair rhs) const -> bool {
  auto [rhsStart, rhsEnd] = rhs;
  return lhs < (rhsEnd - rhsStart);
}

StaticKokkosPool::StaticKokkosPool(size_t numChunks, size_t bytesPerChunks)
    : chunkSize(bytesPerChunks), pool("Memory Pool", numChunks * chunkSize) {
  insertIntoSets({0, numChunks});
}

StaticKokkosPool::StaticKokkosPool(StaticKokkosPool const& other)
    : chunkSize(other.chunkSize), pool("Memory Pool", other.pool.size()) {
  insertIntoSets({0, other.getNumChunks()});
}

auto StaticKokkosPool::insertIntoSets(IndexPair indices)
    -> std::pair<MultiSetBySizeT::iterator, SetByIndexT::iterator> {
  auto setBySizeItr = freeSetBySize.insert(indices);
  auto [setByIndexItr, inserted] = freeSetByIndex.insert(indices);

  assert(inserted);

  return {setBySizeItr, setByIndexItr};
}

void StaticKokkosPool::removeFromSets(IndexPair indices) {
  freeSetBySize.erase(indices);
  freeSetByIndex.erase(indices);
}

auto StaticKokkosPool::allocate(size_t n) -> uint8_t* {
  if (freeSetBySize.empty()) {
    return nullptr;
  }

  // Find the smallest sequence of chunks that can hold numElements
  size_t requestedChunks = getRequiredChunks(n, chunkSize);

  auto freeSetItr = freeSetBySize.lower_bound(requestedChunks);
  if (freeSetItr == freeSetBySize.end()) {
    return nullptr;
  }

  auto [beginIndex, endIndex] = *freeSetItr;

  removeFromSets(*freeSetItr);

  if ((endIndex - beginIndex != requestedChunks)) {
    insertIntoSets({beginIndex + requestedChunks, endIndex});
  }

  uint8_t* ptr = pool.data() + (beginIndex * chunkSize);
  allocations[ptr] = std::make_pair(beginIndex, beginIndex + requestedChunks);

  return ptr;
}

void StaticKokkosPool::deallocate(uint8_t* data) {
  auto allocationsItr = allocations.find(data);
  assert(allocationsItr != allocations.end());
  auto [ptr, chunkIndices] = *allocationsItr;  // [begin, end)

  auto [freeSetBySizeItr, freeSetByIndexItr] = insertIntoSets(chunkIndices);
  allocations.erase(allocationsItr);

  // Merge adjacent free chunks
  if (freeSetByIndexItr != freeSetByIndex.begin()) {
    auto prevItr = std::prev(freeSetByIndexItr);
    auto [prevBeginIndex, prevEndIndex] = *prevItr;

    if (prevEndIndex == chunkIndices.first) {
      removeFromSets(*prevItr);
      removeFromSets(chunkIndices);
      freeSetByIndexItr =
          insertIntoSets({prevBeginIndex, chunkIndices.second}).second;
    }
  }

  if (std::next(freeSetByIndexItr) != freeSetByIndex.end()) {
    auto nextItr = std::next(freeSetByIndexItr);
    auto [nextBeginIndex, nextEndIndex] = *nextItr;
    auto [beginIndex, endIndex] = *freeSetByIndexItr;

    if (chunkIndices.second == nextBeginIndex) {
      removeFromSets(*freeSetByIndexItr);
      removeFromSets(*nextItr);
      insertIntoSets({beginIndex, nextEndIndex});
    }
  }
}

auto StaticKokkosPool::getNumAllocations() const -> unsigned {
  return allocations.size();
}

auto StaticKokkosPool::getNumFreeChunks() const -> unsigned {
  unsigned numFreeChunks = 0;

  for (const auto& [beginIndex, endIndex] : freeSetByIndex) {
    numFreeChunks += endIndex - beginIndex;
  }

  return numFreeChunks;
}

auto StaticKokkosPool::getNumAllocatedChunks() const -> unsigned {
  unsigned numAllocatedChunks = 0;

  for (const auto& [ptr, indices] : allocations) {
    numAllocatedChunks += indices.second - indices.first;
  }

  return numAllocatedChunks;
}

auto StaticKokkosPool::getNumChunks() const -> unsigned {
  return pool.size() / chunkSize;
}

auto StaticKokkosPool::getNumFreeFragments() const -> unsigned {
  return freeSetBySize.size();
}

auto StaticKokkosPool::getRequiredChunks(size_t n, size_t bytesPerChunk)
    -> size_t {
  return (n / bytesPerChunk) + (n % bytesPerChunk ? 1 : 0);
}

void StaticKokkosPool::printDebugInfo() const {
  std::cout
    << "Available Chunks: " << getNumFreeChunks()
    << " Available Fragments: " << getNumFreeFragments()
    << " Allocated Chunks: " << getNumAllocatedChunks()
    << " Allocated Fragments: " << getNumAllocations()
    << " Total Chunks: " << getNumChunks() << std::endl;
}

auto KokkosPool::getChunkSize() const -> size_t { return chunkSize; }

auto KokkosPool::getNumFreeFragments() const -> unsigned {
  unsigned numFreeFragments = 0;

  for (const auto& pool : pools) {
    numFreeFragments += pool.getNumFreeFragments();
  }

  return numFreeFragments;
}

KokkosPool::KokkosPool(size_t initialChunks, size_t bytesPerChunk)
    : chunkSize(bytesPerChunk) {
  pools.emplace_back(initialChunks, chunkSize);
}

KokkosPool::KokkosPool(const KokkosPool& other) : chunkSize(other.chunkSize) {
  for (const auto& pool : other.pools) {
    pools.emplace_back(pool);
  }
}

auto KokkosPool::allocate(size_t n) -> uint8_t* {
  auto current = pools.begin();
  size_t mostAmountOfChunks = 0;

  while (current != pools.end()) {
    uint8_t* ptr = current->allocate(n);
    if (ptr) {
      allocations[ptr] = current;
      return ptr;
    }

    if (current->getNumChunks() > mostAmountOfChunks) {
      mostAmountOfChunks = current->getNumChunks();
    }

    current++;
  }

  size_t requestedChunks = StaticKokkosPool::getRequiredChunks(n, chunkSize);
  size_t amortizedChunkSize = std::max(mostAmountOfChunks * 2, requestedChunks);

  try {
    pools.emplace_back(amortizedChunkSize, chunkSize);
  } catch (const std::runtime_error& e) {

//    for (const auto& pool : pools) {
//      pool.printDebugInfo();
//    }

    std::cout << "Amortization with " << amortizedChunkSize << " chunks failed." << std::endl;
    std::cout << e.what() << std::endl;

    if (amortizedChunkSize <= requestedChunks) {
      throw;
    }

    pools.emplace_back(requestedChunks, chunkSize);
  }

  uint8_t* ptr = pools.back().allocate(n);
  allocations[ptr] = --pools.end();

  return ptr;
}

void KokkosPool::deallocate(uint8_t* data) {
  allocations.at(data)->deallocate(data);
  allocations.erase(data);
}

auto KokkosPool::getNumAllocations() const -> unsigned {
  return allocations.size();
}

auto KokkosPool::getNumFreeChunks() const -> unsigned {
  unsigned numFreeChunks = 0;

  for (const auto& pool : pools) {
    numFreeChunks += pool.getNumFreeChunks();
  }

  return numFreeChunks;
}

auto KokkosPool::getNumAllocatedChunks() const -> unsigned {
  unsigned numAllocatedChunks = 0;

  for (const auto& pool : pools) {
    numAllocatedChunks += pool.getNumAllocatedChunks();
  }

  return numAllocatedChunks;
}

auto KokkosPool::getNumChunks() const -> unsigned {
  unsigned numChunks = 0;

  for (const auto& pool : pools) {
    numChunks += pool.getNumChunks();
  }

  return numChunks;
}

KokkosPool& KokkosPool::operator=(const KokkosPool& other) {
  chunkSize = other.chunkSize;
  pools.clear();
  for (const auto& pool : other.pools) {
    pools.emplace_back(pool);
  }
  return *this;
}

auto KokkosPool::getGlobalPool() -> KokkosPool& {
  if (!s_pool) {
    s_pool = std::make_optional<KokkosPool>(4000, 1'000'000);
  }

  return *s_pool;
}
auto KokkosPool::destroyGlobalPool() -> void { s_pool.reset(); }

}  // namespace Omega_h