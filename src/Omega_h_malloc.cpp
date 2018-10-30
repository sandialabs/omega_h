#include <Omega_h_malloc.hpp>
#include <Omega_h_pool.hpp>
#include <Omega_h_profile.hpp>

namespace Omega_h {

void* device_malloc(std::size_t size) {
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_CUDA
  void* tmp_ptr;
  auto cuda_malloc_size = size;
  if (cuda_malloc_size < 1) cuda_malloc_size = 1;
  auto const err = cudaMalloc(&tmp_ptr, cuda_malloc_size);
  if (err == cudaErrorMemoryAllocation) return nullptr;
  OMEGA_H_CHECK(err == cudaSuccess);
  return tmp_ptr;
#else
  return ::std::malloc(size);
#endif
}

void device_free(void* ptr) {
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_CUDA
  auto const err = cudaFree(ptr);
  OMEGA_H_CHECK(err == cudaSuccess);
#else
  ::std::free(ptr);
#endif
}

void* host_malloc(std::size_t size) {
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_CUDA
  void* tmp_ptr;
  auto cuda_malloc_size = size;
  if (cuda_malloc_size < 1) cuda_malloc_size = 1;
  auto const err = cudaMallocHost(&tmp_ptr, cuda_malloc_size);
  if (err == cudaErrorMemoryAllocation) return nullptr;
  OMEGA_H_CHECK(err == cudaSuccess);
  return tmp_ptr;
#else
  return ::std::malloc(size);
#endif
}

void host_free(void* ptr) {
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_CUDA
  auto const err = cudaFreeHost(ptr);
  OMEGA_H_CHECK(err == cudaSuccess);
#else
  ::std::free(ptr);
#endif
}

static Pool* device_pool = nullptr;
static Pool* host_pool = nullptr;

void enable_pooling() {
#ifdef OMEGA_H_USE_CUDA
  auto const page_size = 2 * 1024 * 1024;
#else
  auto const page_size = 4 * 1024;
#endif
  device_pool = new Pool(page_size, device_malloc, device_free);
  host_pool = new Pool(page_size, host_malloc, host_free);
}

void disable_pooling() {
  delete device_pool;
  delete host_pool;
  device_pool = nullptr;
  host_pool = nullptr;
}

void* maybe_pooled_device_malloc(std::size_t size) {
  if (device_pool) return allocate(*device_pool, size);
  return device_malloc(size);
}

void maybe_pooled_device_free(void* ptr) {
  if (device_pool) deallocate(*device_pool, ptr);
  else device_free(ptr);
}

void* maybe_pooled_host_malloc(std::size_t size) {
  if (host_pool) return allocate(*host_pool, size);
  return host_malloc(size);
}

void maybe_pooled_host_free(void* ptr) {
  if (host_pool) deallocate(*host_pool, ptr);
  else host_free(ptr);
}

}
